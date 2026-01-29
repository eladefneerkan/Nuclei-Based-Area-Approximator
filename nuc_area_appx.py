# Green-to-nucleus assignment

# Outputs per .czi:
#   *_green_assigned_to_nuclei.csv
#   *_GREEN_MASK_CLEAN.tif
#   *_SEEDS_ON_GREEN.tif
#   *_GREEN_ASSIGNED_LABELS.tif
#   *_GREEN_ASSIGNMENT_OUTLINES_ON_GREEN.tif
#   *_GREEN_ASSIGNMENT_OUTLINES_ON_MASK.tif
#   *_UNASSIGNED_GREEN_MASK.tif

from ij import IJ, WindowManager, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.process import ImageConverter, ByteProcessor, ShortProcessor
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
import os, math


# SPEED
BINNING = 2                 # 1 full res (slower), 2 recommended
PRINT_EVERY_GREEN = 200000  # progress log during green assignment

# CHANNELS
GREEN_CHANNEL_INDEX = 1
BLUE_CHANNEL_INDEX  = 2   # nuclei

DO_MAX_Z_PROJECTION = True

# NUCLEI SEGMENTATION (FULL RES)
NUC_BG_SUBTRACT = 30
NUC_BLUR_SIGMA  = 1
NUC_THRESH = "Otsu"
NUC_INVERT_BEFORE_THRESHOLD = True
NUC_FILL_HOLES = True
NUC_DO_WATERSHED = False   # keep OFF 
NUC_MIN_AREA = 10.0
NUC_MAX_AREA = 1e12

# GREEN MASK (BINNED RES)
GRN_BG_SUBTRACT = 50
GRN_BLUR_SIGMA  = 1
GRN_THRESH      = "Li"    # try "Otsu"
GRN_OPEN_CLOSE  = True
GRN_MIN_PARTICLE_AREA = 30.0

# ASSIGNMENT RULES
# Only GREEN pixels get assigned. Empty/background is ignored.
# If a green pixel is farther than this from any nucleus seed, it stays UNASSIGNED.
ASSIGN_MAX_RADIUS_PX = 170   # in BINNED pixels

# Spatial hashing cell size. Use same scale as radius for efficiency.
HASH_CELL_SIZE = ASSIGN_MAX_RADIUS_PX

# OUTPUT
OUT_FOLDER_NAME = "change_me"


# HELPERS

def close_all():
    ids = WindowManager.getIDList()
    if ids:
        for i in ids:
            imp = WindowManager.getImage(i)
            if imp:
                imp.changes = False
                imp.close()

def to_gray8(img):
    try:
        ImageConverter(img).convertToGray8()
    except:
        IJ.run(img, "8-bit", "")
    return img

def open_czi(path):
    IJ.run("Bio-Formats Importer",
           "open=[" + path + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT")
    return WindowManager.getCurrentImage()

def extract_channel(imp, ch_index, title):
    IJ.run(imp, "Duplicate...", "title=%s channel=%d" % (title, ch_index))
    return WindowManager.getCurrentImage()

def max_project(img):
    if DO_MAX_Z_PROJECTION and img.getStackSize() > 1:
        IJ.run(img, "Z Project...", "projection=[Max Intensity]")
        img.changes = False
        img.close()
        return WindowManager.getCurrentImage()
    return img

def ensure_white_objects_on_black(mask_img):
    ip = mask_img.getProcessor()
    mean = ip.getStatistics().mean
    if mean > 127:
        IJ.run(mask_img, "Invert", "")

def rm_get():
    rm = RoiManager.getInstance()
    if rm is None:
        rm = RoiManager()
    return rm

def bin_image(img, binning):
    if binning <= 1:
        return img
    d = img.duplicate()
    IJ.run(d, "Bin...", "x=%d y=%d bin=Average" % (binning, binning))
    return d

def dedup_seeds(seeds):
    seen = set()
    out = []
    for (x, y) in seeds:
        if (x, y) not in seen:
            seen.add((x, y))
            out.append((x, y))
    return out

def draw_seeds_on_image(img, seeds, radius=3):
    d = img.duplicate()
    to_gray8(d)
    ip = d.getProcessor()
    w, h = d.getWidth(), d.getHeight()
    for (x, y) in seeds:
        for yy in range(y - radius, y + radius + 1):
            for xx in range(x - radius, x + radius + 1):
                if 0 <= xx < w and 0 <= yy < h:
                    ip.set(xx, yy, 255)
    d.updateAndDraw()
    return d

def burn_outline_on_image(base_img, outline_bp):
    dup = base_img.duplicate()
    to_gray8(dup)
    ip = dup.getProcessor()
    w = dup.getWidth()
    h = dup.getHeight()
    for y in range(h):
        for x in range(w):
            if outline_bp.get(x, y) == 255:
                ip.set(x, y, 255)
    dup.updateAndDraw()
    return dup
    
def compute_adaptive_radius(seeds, w, h, default_r, min_seeds_for_limit=10, cap=None):
    # If too few seeds, don't clip (avoid circles)
    if len(seeds) < min_seeds_for_limit:
        # big enough to cover the image
        big = int(math.sqrt(w*w + h*h))
        return big if cap is None else min(big, cap)

    return default_r if cap is None else min(default_r, cap)


# Centroids
def nuclei_centroid_seeds(blue_img):
    b = blue_img.duplicate()
    b.setTitle("blue_work")
    to_gray8(b)

    if NUC_BG_SUBTRACT > 0:
        IJ.run(b, "Subtract Background...", "rolling=%d" % NUC_BG_SUBTRACT)
    if NUC_BLUR_SIGMA > 0:
        IJ.run(b, "Gaussian Blur...", "sigma=%d" % NUC_BLUR_SIGMA)
    if NUC_INVERT_BEFORE_THRESHOLD:
        IJ.run(b, "Invert", "")

    IJ.setAutoThreshold(b, NUC_THRESH)
    IJ.run(b, "Convert to Mask", "")
    ensure_white_objects_on_black(b)

    if NUC_FILL_HOLES:
        IJ.run(b, "Fill Holes", "")
    if NUC_DO_WATERSHED:
        IJ.run(b, "Watershed", "")

    rm = rm_get()
    rm.reset()

    IJ.run("Clear Results", "")
    IJ.run("Set Measurements...", "area centroid redirect=None decimal=3")

    IJ.run(b, "Analyze Particles...",
           "size=%s-%s add clear" % (str(NUC_MIN_AREA), str(NUC_MAX_AREA)))

    seeds = []
    w = b.getWidth()
    h = b.getHeight()

    for i in range(rm.getCount()):
        r = rm.getRoi(i)
        st = r.getStatistics()
        x = int(round(st.xCentroid))
        y = int(round(st.yCentroid))
        x = max(0, min(w - 1, x))
        y = max(0, min(h - 1, y))
        seeds.append((x, y))

    print("  nuclei seeds (full-res):", len(seeds))

    b.changes = False
    b.close()
    return seeds


# GREEN MASK (BINNED) (binary)
def green_mask_clean_at_binned_res(green_b):
    g = green_b.duplicate()
    g.setTitle("green_work_binned")
    to_gray8(g)

    if GRN_BG_SUBTRACT > 0:
        IJ.run(g, "Subtract Background...", "rolling=%d" % GRN_BG_SUBTRACT)
    if GRN_BLUR_SIGMA > 0:
        IJ.run(g, "Gaussian Blur...", "sigma=%d" % GRN_BLUR_SIGMA)

    IJ.setAutoThreshold(g, GRN_THRESH)
    IJ.run(g, "Convert to Mask", "")
    ensure_white_objects_on_black(g)

    if GRN_OPEN_CLOSE:
        IJ.run(g, "Open", "")
        IJ.run(g, "Close", "")

    IJ.run("Clear Results", "")
    IJ.run("Set Measurements...", "area redirect=None decimal=3")
    IJ.run(g, "Analyze Particles...",
           "size=%s-Infinity show=Masks clear" % str(GRN_MIN_PARTICLE_AREA))

    g.changes = False
    g.close()
    cleaned = WindowManager.getCurrentImage()
    cleaned.setTitle("GREEN_MASK_CLEAN")
    return cleaned


# Spatial hash for seeds (fast nearest within radius)
def build_seed_hash(seeds, cell_size):
    # map (cx,cy) -> list of (sx,sy,label)
    H = {}
    for i, (sx, sy) in enumerate(seeds):
        label = i + 1
        cx = int(sx // cell_size)
        cy = int(sy // cell_size)
        key = (cx, cy)
        if key not in H:
            H[key] = []
        H[key].append((sx, sy, label))
    return H

def nearest_seed_within_radius(x, y, H, cell_size, max_r):
    if max_r <= 0:
        return 0

    cx = int(x // cell_size)
    cy = int(y // cell_size)
    r2 = max_r * max_r

    best_lab = 0
    best_d2 = r2 + 1

    # check this cell and neighbors (3x3)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            key = (cx + dx, cy + dy)
            if key not in H:
                continue
            for (sx, sy, lab) in H[key]:
                ddx = x - sx
                ddy = y - sy
                d2 = ddx * ddx + ddy * ddy
                if d2 <= r2 and d2 < best_d2:
                    best_d2 = d2
                    best_lab = lab

    return best_lab

# Assign GREEN pixels
def assign_green_to_nuclei(gmask_b, seeds_b, max_r, cell_size):
    ip = gmask_b.getProcessor()
    w, h = gmask_b.getWidth(), gmask_b.getHeight()
    n = w * h

    # Build hash
    H = build_seed_hash(seeds_b, cell_size)

    labels = [0] * n
    green_total = 0
    unassigned_green = 0

    # Per nucleus counts
    counts = [0] * (len(seeds_b) + 1)  # 1..N

    visited = 0
    for y in range(h):
        row = y * w
        for x in range(w):
            if ip.get(x, y) <= 0:
                continue

            green_total += 1
            lab = nearest_seed_within_radius(x, y, H, cell_size, max_r)
            if lab <= 0:
                unassigned_green += 1
                continue

            labels[row + x] = lab
            counts[lab] += 1

            visited += 1
            if PRINT_EVERY_GREEN > 0 and (visited % PRINT_EVERY_GREEN == 0):
                print("  green assignment progress:", visited, "assigned pixels")

    print("  green_total_px:", green_total, "unassigned_green_px:", unassigned_green)
    return labels, counts, green_total, unassigned_green


# Outline boundaries of assignment labels (only where label>0)
def outline_from_assignment_labels(width, height, labels):
    out = ByteProcessor(width, height)
    for y in range(height):
        row = y * width
        for x in range(width):
            idx = row + x
            v = labels[idx]
            if v == 0:
                out.set(x, y, 0)
                continue
            boundary = False
            if x > 0 and labels[idx - 1] != v: boundary = True
            elif x < width - 1 and labels[idx + 1] != v: boundary = True
            elif y > 0 and labels[idx - width] != v: boundary = True
            elif y < height - 1 and labels[idx + width] != v: boundary = True
            out.set(x, y, 255 if boundary else 0)
    return out

def labels_to_short_image(width, height, labels, title):
    sp = ShortProcessor(width, height)
    spix = sp.getPixels()  # short[]
    for i in range(width * height):
        spix[i] = labels[i]
    return ImagePlus(title, sp)

def unassigned_green_mask(gmask_b, assigned_labels):
    ip = gmask_b.getProcessor()
    w, h = gmask_b.getWidth(), gmask_b.getHeight()
    out = ByteProcessor(w, h)
    for y in range(h):
        row = y * w
        for x in range(w):
            if ip.get(x, y) > 0 and assigned_labels[row + x] == 0:
                out.set(x, y, 255)
            else:
                out.set(x, y, 0)
    return ImagePlus("UNASSIGNED_GREEN_MASK", out)


# main
dc = DirectoryChooser("Choose folder with CZI files")
inputDir = dc.getDirectory()
if inputDir is None:
    raise SystemExit("No folder selected")

outputDir = os.path.join(inputDir, OUT_FOLDER_NAME)
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

files = [f for f in os.listdir(inputDir) if f.lower().endswith(".czi")]
if not files:
    raise SystemExit("No .czi files found in: " + inputDir)

for f in files:
    close_all()
    print("Processing:", f)

    imp = open_czi(os.path.join(inputDir, f))
    green = extract_channel(imp, GREEN_CHANNEL_INDEX, "green")
    blue  = extract_channel(imp, BLUE_CHANNEL_INDEX,  "blue")

    green = max_project(green)
    blue  = max_project(blue)

    # seeds full-res
    seeds_full = nuclei_centroid_seeds(blue)
    if len(seeds_full) == 0:
        print("  WARNING: no nuclei found; skipping", f)
        continue

    # bin green for speed, and scale seeds to binned coordinates
    green_b = bin_image(green, BINNING)
    w = green_b.getWidth()
    h = green_b.getHeight()

    if BINNING > 1:
        seeds_b = [(int(round(x / float(BINNING))), int(round(y / float(BINNING)))) for (x, y) in seeds_full]
    else:
        seeds_b = list(seeds_full)

    seeds_b = dedup_seeds(seeds_b)
    print("  seeds after bin+dedup:", len(seeds_b))

    # clean green mask at binned res
    gmask_b = green_mask_clean_at_binned_res(green_b)

    base = os.path.splitext(f)[0]

    # QC: seeds on green
    seed_qc = draw_seeds_on_image(green_b, seeds_b, radius=3)
    seed_qc.setTitle("SEEDS_ON_GREEN")
    FileSaver(seed_qc).saveAsTiff(os.path.join(outputDir, base + "_SEEDS_ON_GREEN.tif"))
    seed_qc.changes = False
    seed_qc.close()

    # Assign green pixels to nearest nucleus
    # Adaptive radius to avoid artificial circular cutoffs :)
    R = compute_adaptive_radius(
        seeds_b,
        w,
        h,
        ASSIGN_MAX_RADIUS_PX,
        min_seeds_for_limit=10
    )

    labels, counts, green_total, unassigned_green = assign_green_to_nuclei(
        gmask_b,
        seeds_b,
        R,
        HASH_CELL_SIZE
    )

    print("  using radius:", R)



    # Label image (only green pixels have labels)
    lab_imp = labels_to_short_image(w, h, labels, "GREEN_ASSIGNED_LABELS")

    # Outline boundaries of green assignment
    outline_bp = outline_from_assignment_labels(w, h, labels)
    outlines_on_green = burn_outline_on_image(green_b, outline_bp)
    outlines_on_green.setTitle("GREEN_ASSIGNMENT_OUTLINES_ON_GREEN")
    outlines_on_mask  = burn_outline_on_image(gmask_b, outline_bp)
    outlines_on_mask.setTitle("GREEN_ASSIGNMENT_OUTLINES_ON_MASK")

    # Unassigned green QC
    unass = unassigned_green_mask(gmask_b, labels)
    unass.setTitle("UNASSIGNED_GREEN_MASK")

    # CSV per nucleus
    table = ResultsTable()
    for i, (sx, sy) in enumerate(seeds_b):
        nid = i + 1
        table.incrementCounter()
        table.addValue("Nucleus_ID", nid)
        table.addValue("SeedX", sx)
        table.addValue("SeedY", sy)
        table.addValue("GreenArea_px_assigned", counts[nid])
        frac_total_green = (100.0 * counts[nid] / green_total) if green_total > 0 else 0.0
        table.addValue("Fraction_of_total_green_%", frac_total_green)

    # also add a summary row for unassigned
    table.incrementCounter()
    table.addValue("Nucleus_ID", 0)
    table.addValue("SeedX", -1)
    table.addValue("SeedY", -1)
    table.addValue("GreenArea_px_assigned", unassigned_green)
    frac_un = (100.0 * unassigned_green / green_total) if green_total > 0 else 0.0
    table.addValue("Fraction_of_total_green_%", frac_un)

    csv_path = os.path.join(outputDir, base + "_green_assigned_to_nuclei.csv")
    table.save(csv_path)

    # Save outputs
    FileSaver(gmask_b).saveAsTiff(os.path.join(outputDir, base + "_GREEN_MASK_CLEAN.tif"))
    FileSaver(lab_imp).saveAsTiff(os.path.join(outputDir, base + "_GREEN_ASSIGNED_LABELS.tif"))
    FileSaver(outlines_on_green).saveAsTiff(os.path.join(outputDir, base + "_GREEN_ASSIGNMENT_OUTLINES_ON_GREEN.tif"))
    FileSaver(outlines_on_mask).saveAsTiff(os.path.join(outputDir, base + "_GREEN_ASSIGNMENT_OUTLINES_ON_MASK.tif"))
    FileSaver(unass).saveAsTiff(os.path.join(outputDir, base + "_UNASSIGNED_GREEN_MASK.tif"))

    print("  Saved:", csv_path)

    # cleanup
    for img in [imp, green, blue, green_b, gmask_b, lab_imp, outlines_on_green, outlines_on_mask, unass]:
        try:
            img.changes = False
            img.close()
        except:
            pass

print("DONE :)")
