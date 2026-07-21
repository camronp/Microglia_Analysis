// added border margin of 4 pixels
// added cropped image naming

var threshold_algorithm, min_cell_size, max_cell_size;
var subtract_bg_radius, smooth_radius, close_gaps, fill_detection_holes;
var split_touching, split_dynamic, mask_grow, padding, fallback_um_per_px;
var input_folder, output_folder, crops_folder, csv_folder, maps_folder;
var images, ni, grand_total_cells, preset, accept_all;
var accept_all_batch, border_margin;
var accept_all_batch, border_margin, saved_skeleton_example;

Dialog.create("Microglia Segmentation");
Dialog.addMessage("Pick a starting preset. You can preview and change it\nfor each image before it's processed.");
Dialog.addChoice("Starting image type",
    newArray("Standard", "Dim / noisy", "Small cells", "ImageXpress"), "ImageXpress");
Dialog.addString("Only process files containing (e.g. w2; blank = all)", "w2");
Dialog.show();
preset = Dialog.getChoice();
channel_filter = Dialog.getString();


fill_detection_holes = true;
close_gaps = 4;           
split_touching = false;   
split_dynamic  = 2;
accept_all = false;
accept_all_batch = false;
border_margin = 4;   // px; a cell within this many px of the edge counts as "touching". Raise if edge cells still slip through.
saved_skeleton_example = false;

applyPreset(preset);      


if (preset == "Small cells")            fallback_um_per_px = 4 * (184.65 / 512);
else if (preset == "ImageXpress") fallback_um_per_px = 0.65;  
else                                        fallback_um_per_px = 184.65 / 512;


input_folder  = getDirectory("Choose the INPUT folder (original images)");
output_folder = getDirectory("Choose the OUTPUT folder (results go here)");

crops_folder = output_folder + "output_images/";
csv_folder   = output_folder + "analysis_csv/";
maps_folder  = output_folder + "cell_maps/";
if (!File.exists(crops_folder)) File.makeDirectory(crops_folder);
if (!File.exists(csv_folder))   File.makeDirectory(csv_folder);
if (!File.exists(maps_folder))  File.makeDirectory(maps_folder);


file_list = getFileList(input_folder);
images = newArray();
ni = 0;
channel_lc = toLowerCase(channel_filter);
for (fi = 0; fi < file_list.length; fi++) {
    lname = toLowerCase(file_list[fi]);
    is_tif = endsWith(lname, ".tif") || endsWith(lname, ".tiff");
    matches_channel = (channel_lc == "") || (indexOf(lname, channel_lc) >= 0);
    if (is_tif && matches_channel) {
        images[ni] = file_list[fi];
        ni++;
    }
}

print("\\Clear");
print("=== Microglia Segment + Analysis ===");
print("Input : " + input_folder);
print("Output: " + output_folder);
if (channel_filter != "") print("Channel filter: '" + channel_filter + "'");
print("Found " + ni + " image(s) to process");
print("");

if (ni == 0) {
    showMessage("No images", "No .tif images matching channel '" + channel_filter +
        "' found in:\n" + input_folder + "\n\n(Clear the channel field to process all images.)");
    exit;
}


grand_total_cells = 0;
for (img = 0; img < ni; img++) {
    print("=== IMAGE " + (img+1) + "/" + ni + ": " + images[img] + " ===");
    processImage(images[img]);
    print("");
}

print("=== ALL DONE ===");
print("Total cells analyzed: " + grand_total_cells);
print("Crops: " + crops_folder);
print("Maps : " + maps_folder);
print("CSVs : " + csv_folder);
showMessage("Complete",
    "Processed " + ni + " image(s), " + grand_total_cells + " cells total.\n\n" +
    "Crops: output_images/ (<image_id>_cell<N>.tif)\n" +
    "Maps:  cell_maps/     (<image_id>_map.tif, numbered)\n" +
    "CSVs:  analysis_csv/ (R-pipeline ready)");




function processImage(image_filename) {
    open(input_folder + image_filename);
    original_title = getTitle();
    getDimensions(width, height, channels, slices, frames);
    print("  " + width + " x " + height + " px");

    accept_all = tuneDetection(original_title, width, height, image_filename);

    detected = roiManager("count");
    selectWindow(original_title);
    Overlay.remove;
    print("  Detected: " + detected + " candidate cells");

    if (detected == 0) {
        print("  (skipping image - no cells)");
        close(original_title);
        return;
    }

 
    approved_roi = newArray();
    approved_x   = newArray();
    approved_y   = newArray();
    approved_w   = newArray();
    approved_h   = newArray();
    napproved = 0;

    if (accept_all) {
   
        for (i = 0; i < detected; i++) {
            roiManager("select", i);
            getSelectionBounds(x, y, w, h);
            approved_roi[napproved] = i;
            approved_x[napproved] = x;
            approved_y[napproved] = y;
            approved_w[napproved] = w;
            approved_h[napproved] = h;
            napproved++;
        }
    } else {
        total = detected;
        i = 0;
        while (i < total) {
            selectWindow(original_title);
            roiManager("select", i);
            getSelectionBounds(x, y, w, h);
            Overlay.remove;
            Overlay.addSelection("yellow");
            Overlay.show;

            Dialog.create("Review cell " + (i+1) + " of " + total);
            Dialog.addMessage(images[img]);
            Dialog.addMessage("Position X=" + x + ", Y=" + y + "   Size " + w + "x" + h);
            Dialog.addRadioButtonGroup("Action", newArray("KEEP", "SKIP", "SPLIT"), 1, 1, "KEEP");
            Dialog.show();
            act = Dialog.getRadioButton();
            Overlay.remove;

            if (act == "KEEP") {
                approved_roi[napproved] = i;
                approved_x[napproved] = x;
                approved_y[napproved] = y;
                approved_w[napproved] = w;
                approved_h[napproved] = h;
                napproved++;
            } else if (act == "SPLIT") {
                splitCell(i, original_title, width, height);
                total = roiManager("count");
            }
            i++;
        }
    }
    print("  Approved " + napproved + " cells");

    if (napproved == 0) {
        close(original_title);
        return;
    }

   
    getPixelSize(cal_unit, pw, ph);
    if (pw == 0 || cal_unit == "pixel") {
        selectWindow(original_title);
        setVoxelSize(fallback_um_per_px, fallback_um_per_px, 1, "um");
        getPixelSize(cal_unit, pw, ph);
    }

  
    getDateAndTime(yr, mo, dow, dom, hr, mn, sc, ms);
    date_str = "" + yr + (mo+1) + dom;
    time_str = "" + hr + mn + sc;
    base_name  = substring(image_filename, 0, lengthOf(image_filename) - 4);
    clean_base = replace(base_name, "[^A-Za-z0-9]", "_");
    csv_path = csv_folder + "morphology_results_" + clean_base + "_" + date_str + time_str + ".csv";

    header = "cell_name,area_um2,convexhull_area_um2,perimeter_um,convexhull_perimeter_um," +
             "ramification_index,circularity,solidity,n_branches,n_junctions,tree_length_um," +
             "avg_branch_length_um,max_branch_length_um,n_tips,n_triple_points,n_quadruple_points," +
             "aspect_ratio,polarity_offset_um,centroid_x_um,centroid_y_um";
    File.saveString(header + "\n", csv_path);

    setBatchMode(true);


    selectWindow(original_title);
    newImage("cell_map", "8-bit black", width, height, 1);
    setColor(255);
    for (m = 0; m < napproved; m++) {
        roiManager("select", approved_roi[m]);
        fill();
    }
    run("Select None");

    Overlay.remove;
    for (m = 0; m < napproved; m++) {
        roiManager("select", approved_roi[m]);
        Overlay.addSelection("cyan");
    }
    run("Select None");

    setFont("SansSerif", 30, "bold");
    setColor("gray");
    for (m = 0; m < napproved; m++) {
        lx = approved_x[m] + approved_w[m] / 2 - 8;
        ly = approved_y[m] + approved_h[m] / 2 + 8;
        Overlay.drawString("" + (m + 1), lx, ly);
    }
    Overlay.show;
    run("Flatten");
    map_title = getTitle();
    saveAs("TIFF", maps_folder + base_name + "_map.tif");
    close(map_title);
    close("cell_map");
    selectWindow(original_title);

    for (k = 0; k < napproved; k++) {
        x = approved_x[k];  y = approved_y[k];
        w = approved_w[k];  h = approved_h[k];
        ridx = approved_roi[k];

     
        selectWindow(original_title);
        roiManager("deselect");
        run("Select None");
        run("Duplicate...", "title=cell_work");
        roiManager("select", ridx);
        run("Enlarge...", "enlarge=" + mask_grow);
        setBackgroundColor(0, 0, 0);
        run("Clear Outside");

        run("Select None");
        cx = maxOf(0, x - padding);
        cy = maxOf(0, y - padding);
        cw = minOf(w + 2*padding, width  - cx);
        ch = minOf(h + 2*padding, height - cy);
        makeRectangle(cx, cy, cw, ch);
        run("Crop");

        crop_name = base_name + "_cell" + (k + 1) + ".tif";
        saveAs("TIFF", crops_folder + crop_name);
        crop_title = getTitle();


        if (bitDepth() != 8) run("8-bit");
        if (smooth_radius > 0) run("Median...", "radius=" + smooth_radius);
        setAutoThreshold(threshold_algorithm + " dark");
        run("Create Selection");
        resetThreshold();
        run("Set Measurements...",
            "area perimeter shape centroid center redirect=None decimal=4");
        run("Clear Results");
        run("Measure");
        area     = getResult("Area", 0);
        perim    = getResult("Perim.", 0);
        circ     = getResult("Circ.", 0);
        ar       = getResult("AR", 0);
        solidity = getResult("Solidity", 0);
        gx = getResult("X", 0);   gy = getResult("Y", 0);
        mx = getResult("XM", 0);  my = getResult("YM", 0);
        polarity = sqrt((gx-mx)*(gx-mx) + (gy-my)*(gy-my));

       
        run("Restore Selection");
        run("Convex Hull");
        run("Measure");
        convex_area  = getResult("Area", 1);
        convex_perim = getResult("Perim.", 1);
        if (convex_perim > 0) ramification = perim / convex_perim; else ramification = 0;

        selectWindow(crop_title);
        run("Select None");
        run("Duplicate...", "title=skel_work");
        if (bitDepth() != 8) run("8-bit");
        setOption("BlackBackground", true);
        setAutoThreshold(threshold_algorithm + " dark");
        run("Convert to Mask");
        run("Fill Holes");
        run("Skeletonize");
        run("Clear Results");
        run("Analyze Skeleton (2D/3D)", "prune=none");

        n_branches=0; n_junctions=0; n_tips=0; n_triple=0; n_quad=0;
        avg_branch=0; max_branch=0; tree_len=0;
        if (nResults > 0) {
            best=-1; bestBr=-1;
            for (r=0; r<nResults; r++) {
                br = getResult("# Branches", r);
                if (!isNaN(br) && br > bestBr) { bestBr = br; best = r; }
            }
            if (best >= 0) {
                n_branches  = getResult("# Branches", best);
                n_junctions = getResult("# Junctions", best);
                n_tips      = getResult("# End-point voxels", best);
                n_triple    = getResult("# Triple points", best);
                n_quad      = getResult("# Quadruple points", best);
                avg_branch  = getResult("Average Branch Length", best);
                max_branch  = getResult("Maximum Branch Length", best);
                if (!isNaN(n_branches) && !isNaN(avg_branch)) tree_len = n_branches * avg_branch;
            }
        }
        close("skel_work");
        if (!saved_skeleton_example && isOpen("Tagged skeleton")) {
    selectWindow("Tagged skeleton");
    saveAs("TIFF", output_folder + "example_tagged_skeleton_" + base_name + "_cell" + (k + 1) + ".tif");
    saved_skeleton_example = true;
}
if (isOpen("Tagged skeleton")) close("Tagged skeleton");
        if (isOpen("Tagged skeleton")) close("Tagged skeleton");
        if (isOpen("Longest shortest paths")) close("Longest shortest paths");

        
        centroid_x = cx * pw + gx;
        centroid_y = cy * ph + gy;

        row = "\"" + crop_name + "\"," +
              d2s(area,4) + "," + d2s(convex_area,4) + "," +
              d2s(perim,4) + "," + d2s(convex_perim,4) + "," +
              d2s(ramification,4) + "," + d2s(circ,4) + "," + d2s(solidity,4) + "," +
              n_branches + "," + n_junctions + "," + d2s(tree_len,4) + "," +
              d2s(avg_branch,4) + "," + d2s(max_branch,4) + "," +
              n_tips + "," + n_triple + "," + n_quad + "," +
              d2s(ar,4) + "," + d2s(polarity,4) + "," +
              d2s(centroid_x,4) + "," + d2s(centroid_y,4);
        File.append(row, csv_path);

        
        titles = getList("image.titles");
        for (t = 0; t < titles.length; t++) {
            if (titles[t] != original_title) { selectWindow(titles[t]); close(); }
        }
    }
    setBatchMode(false);

    roiManager("reset");
    if (isOpen(original_title)) { selectWindow(original_title); close(); }
    grand_total_cells += napproved;
    print("  CSV: " + csv_path);
}


function applyPreset(name) {
    if (name == "Standard 40x") {
        threshold_algorithm = "Li";
        min_cell_size = 150;   max_cell_size = 5000;
        subtract_bg_radius = 50; smooth_radius = 1;
        mask_grow = 4;         padding = 30;
    } else if (name == "Dim / noisy") {
        threshold_algorithm = "Yen";
        min_cell_size = 300;   max_cell_size = 6000;
        subtract_bg_radius = 80; smooth_radius = 1;
        mask_grow = 6;         padding = 30;
    } else if (name == "Small cells") {
        threshold_algorithm = "Triangle";
        min_cell_size = 30;    max_cell_size = 1200;
        subtract_bg_radius = 40; smooth_radius = 1;
        mask_grow = 2;         padding = 12;
    } else if (name == "ImageXpress") {
        
        threshold_algorithm = "Yen";
        min_cell_size = 2000;  max_cell_size = 200000;  
        subtract_bg_radius = 0;  smooth_radius = 4;      
        mask_grow = 8;         padding = 60;
    }
}


function detectCells(original_title, width, height) {
    selectWindow(original_title);
    run("Select None");
    run("Duplicate...", "title=detect_mask");
    if (bitDepth() != 8) run("8-bit");
    if (subtract_bg_radius > 0) run("Subtract Background...", "rolling=" + subtract_bg_radius);
    if (smooth_radius > 0)      run("Median...", "radius=" + smooth_radius);
    setAutoThreshold(threshold_algorithm + " dark");
    run("Convert to Mask");
    if (close_gaps > 0) {
        setOption("BlackBackground", true);
        for (cg = 0; cg < close_gaps; cg++) run("Dilate");
        for (cg = 0; cg < close_gaps; cg++) run("Erode");
    }
    if (fill_detection_holes) run("Fill Holes");

    roiManager("reset");
    if (split_touching) {
        selectWindow("detect_mask");
        run("Distance Transform Watershed",
            "distances=[Borgefors (3,4)] output=[16 bits] normalize" +
            " dynamic=" + split_dynamic + " connectivity=4");
        ws_title = getTitle();
        run("Remove Border Labels", "left right top bottom");
		ws_title = getTitle();
        run("Conversions...", "scale");
        setMinAndMax(0, 255);
        run("8-bit");
        getHistogram(values, counts, 256);
        for (L = 1; L <= 255; L++) {
            if (counts[L] >= min_cell_size && counts[L] <= max_cell_size) {
                selectWindow(ws_title);
                setThreshold(L, L);
                run("Create Selection");
                resetThreshold();
                roiManager("add");
            }
        }
        close(ws_title);
    } else {
        run("Analyze Particles...", "size=" + min_cell_size + "-" + max_cell_size + " exclude add");
    }
    removeBorderROIs(width, height);

    close("detect_mask");
    selectWindow(original_title);
}


function removeBorderROIs(width, height) {

    if (isOpen("detect_mask")) selectWindow("detect_mask");

    n = roiManager("count");
    removed = 0;
    for (r = n - 1; r >= 0; r--) {
        roiManager("select", r);
        getSelectionBounds(bx, by, bw, bh);
        touches_left   = (bx <= border_margin);
        touches_top    = (by <= border_margin);
        touches_right  = ((bx + bw) >= (width  - border_margin));
        touches_bottom = ((by + bh) >= (height - border_margin));
        if (touches_left || touches_top || touches_right || touches_bottom) {
            roiManager("delete");
            removed++;
        }
    }
    roiManager("deselect");
    print("    Border check (image " + width + "x" + height + ", margin=" + border_margin +
          "px): " + n + " candidate(s) -> " + removed + " removed, " + (n - removed) + " kept");
}


function tuneDetection(original_title, width, height, image_filename) {
    if (accept_all_batch) { detectCells(original_title, width, height); return true; }
    accepted = false;
    accept_all = false;
    while (!accepted) {
       
        detectCells(original_title, width, height);
        n = roiManager("count");
        selectWindow(original_title);
        Overlay.remove;
        for (r = 0; r < n; r++) { roiManager("select", r); Overlay.addSelection("cyan"); }
        Overlay.show;
        run("Select None");

        Dialog.create("Preview: " + n + " cells detected");
        Dialog.addMessage(image_filename);
        Dialog.addMessage("Method '" + threshold_algorithm + "' found " + n + " cells (cyan outlines).");
        Dialog.addMessage("Switch preset or open Advanced, then preview again.\n" +
                          "When the outlines look right, choose how to proceed:");
        Dialog.addChoice("Preset",
            newArray("(keep current)", "Standard 40x", "Dim / noisy 40x", "Small cells 10x", "ImageXpress"),
            "(keep current)");
        Dialog.addCheckbox("Advanced settings", false);
        Dialog.addCheckbox("Accept - review each cell (KEEP/SKIP/SPLIT)", false);
        Dialog.addCheckbox("Accept ALL - keep every cell, skip review", false);
        Dialog.addCheckbox("Accept ALL for ALL remaining images", false);
        Dialog.show();

        chosen   = Dialog.getChoice();
        advanced = Dialog.getCheckbox();
        accept_review = Dialog.getCheckbox();
        accept_all    = Dialog.getCheckbox();
        accept_all_batch = Dialog.getCheckbox();
if (accept_all_batch) { accept_all = true; accepted = true; }
        accepted = accept_review || accept_all;

     
        if (chosen != "(keep current)") { applyPreset(chosen); accepted = false; accept_all = false; }
        if (advanced) { advancedDialog(); accepted = false; accept_all = false; }
    }
    Overlay.remove;
    return accept_all;
}

function advancedDialog() {
    Dialog.create("Advanced settings");
    Dialog.addChoice("Threshold method",
        newArray("Li", "Yen", "Triangle", "RenyiEntropy", "MaxEntropy", "Otsu", "Default", "Moments"),
        threshold_algorithm);
    Dialog.addNumber("Min cell size (px)", min_cell_size);
    Dialog.addNumber("Max cell size (px)", max_cell_size);
    Dialog.addNumber("Background subtract radius (px, 0=off)", subtract_bg_radius);
    Dialog.addNumber("Smooth / median radius (px, 0=off)", smooth_radius);
    Dialog.addNumber("Gap-bridging (close, 0=off)", close_gaps);
    Dialog.addNumber("Border exclusion margin (px, 0=must literally touch edge)", border_margin);
    Dialog.show();
    threshold_algorithm = Dialog.getChoice();
    min_cell_size      = Dialog.getNumber();
    max_cell_size      = Dialog.getNumber();
    subtract_bg_radius = Dialog.getNumber();
    smooth_radius      = Dialog.getNumber();
    close_gaps         = Dialog.getNumber();
    border_margin      = Dialog.getNumber();
}


function splitCell(idx, original_title, width, height) {
    selectWindow(original_title);
    setTool("line");
    roiManager("select", idx);
    waitForUser("Split cell",
        "Draw a straight line across where the cells meet, then click OK.\n" +
        "(Line tool is already selected.)");
    if (selectionType() != 5) {   // 5 = straight line
        showMessage("No line", "No straight line drawn cell left unsplit.");
        return;
    }
    getLine(x1, y1, x2, y2, lw);

    newImage("splitmask", "8-bit black", width, height, 1);
    roiManager("select", idx);
    setColor(255);
    fill();
    run("Select None");

    setColor(0);
    setLineWidth(3);
    drawLine(x1, y1, x2, y2);
    setThreshold(128, 255);
    run("Analyze Particles...", "size=" + min_cell_size + "-Infinity pixel add");
    resetThreshold();
    close("splitmask");
    selectWindow(original_title);
}
