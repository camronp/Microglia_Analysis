// added border margin of 4 pixels
// added cropped image naming
// v1.5: also exports the per-cell binary segmentation mask (same threshold
//       used for morphology measurements) + a mask_index CSV linking each
//       mask back to its source image + crop rectangle. This lets the
//       separate "Channel Intensity Quantification" pipeline apply the
//       exact same cell mask to other channels (w1, w3, w4, etc.).
// v1.6: replaced the single "first tagged skeleton" save with a random
//       sample of N tagged-skeleton example images (user-specified count),
//       chosen via reservoir sampling across the whole batch, saved into
//       their own tagged_skeleton_examples/ subfolder.
// v1.7: the saved skeleton examples now show the tagged (color-coded)
//       skeleton composited directly on top of that cell's cropped image,
//       instead of the skeleton alone on a black background.
// 		 output folder is now auto-detected from the macro's own file
//       location (no more "choose output folder" prompt); and all
//       per-image morphology CSVs are now also merged into a single
//       merged_morphology_data.csv (in analysis_csv/, alongside the
//       per-image CSVs), with the same derived columns (sample_id,
//       image_id, cell_number, project, well_position, site,
//       wavelength, image, ramification_index_2d) that data_helpers.R's
//       load_morphology_data() adds during its own post-processing merge.
// 		 added a fixed_output_folder override (skips auto-detection
//       entirely) and a sanity check on the auto-detected macro path,
//       since Fiji's Script Editor "Run" button doesn't always report
//       getInfo("macro.filepath") correctly (can return an unrelated,
//       stale path instead of this macro's own file).
// 		 the merged CSV now always saves directly next to this macro
//       file (not inside analysis_csv/), and its filename now includes
//       the input folder's name (as an experiment tag) and the date/
//       time the run started, e.g. merged_morphology_data_<experiment>_
//       <date><time>.csv.
// 		 added a dedicated fixed_macro_folder override (separate from
//       fixed_output_folder) for the merged CSV specifically, since
//       getInfo("macro.filepath") auto-detection isn't reliable when run
//       from Fiji's Script Editor; the log now also states outright
//       whether it auto-detected the macro's folder or fell back to
//       output_folder, so a wrong location is easy to spot and fix.

var threshold_algorithm, min_cell_size, max_cell_size;
var subtract_bg_radius, smooth_radius, close_gaps, fill_detection_holes;
var split_touching, split_dynamic, mask_grow, padding, fallback_um_per_px;
var input_folder, output_folder, crops_folder, csv_folder, maps_folder, masks_folder;
var images, ni, grand_total_cells, preset, accept_all;
var accept_all_batch, border_margin;
var accept_all_batch, border_margin;
var num_skeleton_examples, skeleton_examples_folder, skeleton_seen_count, reservoir_files;
var merged_csv_path;

Dialog.create("Microglia Segmentation");
Dialog.addMessage("Pick a starting preset. You can preview and change it\nfor each image before it's processed.");
Dialog.addChoice("Starting image type",
    newArray("Standard", "Dim / noisy", "Small cells", "ImageXpress"), "ImageXpress");
Dialog.addString("Only process files containing (e.g. w2; blank = all)", "w2");
Dialog.addNumber("Random tagged-skeleton examples to save (0 = none)", 5);
Dialog.show();
preset = Dialog.getChoice();
channel_filter = Dialog.getString();
num_skeleton_examples = Dialog.getNumber();


fill_detection_holes = true;
close_gaps = 4;           
split_touching = false;   
split_dynamic  = 2;
accept_all = false;
accept_all_batch = false;
border_margin = 4;   // raise to 5 if it is still thresholding cells on the edge
skeleton_seen_count = 0;

// If you always run this from the same project folder, set this to that
// path once (e.g. "C:/Cam/Local_R_work/microglia_analysis/") and every run
// will use it with zero prompts and zero guessing. Leave it "" to fall back
// to auto-detecting the macro's own folder, and then to a folder picker if
// that fails too (see the fixed_output_folder block below for why the
// auto-detect step can't always be trusted).
fixed_output_folder = "";

// Same idea, but specifically for where the merged CSV gets saved (see
// merged_csv_path below) - kept separate from fixed_output_folder because
// you might want the big per-image outputs on a data drive but the merged
// CSV always sitting right next to this script. getInfo("macro.filepath")
// does not reliably report this macro's own path when run from Fiji's
// Script Editor (a known Fiji limitation, not fixable from inside the
// macro) - if the log says it fell back to output_folder, set this to the
// folder this .ijm file actually lives in and it'll be used directly, no
// detection involved.
fixed_macro_folder = "C:/Cam/Local_R_work/microglia_analysis/";

applyPreset(preset);      


if (preset == "Small cells")            fallback_um_per_px = 4 * (184.65 / 512);
else if (preset == "ImageXpress") fallback_um_per_px = 0.65;  
else                                        fallback_um_per_px = 184.65 / 512;


input_folder  = getDirectory("Choose the input folder (original images)");

// Always try to detect the macro's own folder - used below as a fallback
// for output_folder, and separately as the fixed home for the merged CSV
// (see merged_csv_path further down), regardless of where output_folder
// itself ends up (e.g. if fixed_output_folder points somewhere else).
macro_path = getInfo("macro.filepath");
looks_valid = (macro_path != "") && endsWith(toLowerCase(macro_path), ".ijm");
if (looks_valid) macro_folder = File.getDirectory(macro_path);

if (fixed_output_folder != "") {
    output_folder = fixed_output_folder;
    if (!endsWith(output_folder, "/") && !endsWith(output_folder, "\\")) output_folder = output_folder + "/";
} else if (looks_valid) {
    output_folder = macro_folder;
} else {
    // Neither a fixed folder nor a reliably auto-detected macro path
    // (see the note on Fiji's Script Editor above fixed_output_folder) -
    // fall back to a folder picker instead of writing output to the wrong
    // place. If you hit this, set fixed_output_folder above instead.
    output_folder = getDirectory("Choose the output folder (results go here)");
}
if (!File.exists(output_folder)) File.makeDirectory(output_folder);

// Where the merged CSV goes: fixed_macro_folder always wins if set (most
// reliable - use it if the log below ever says auto-detect failed). Then
// the auto-detected macro path, if it looked valid. Otherwise falls back
// to output_folder, with a console line explaining why, so this is
// diagnosable from the log instead of just silently landing somewhere else.
if (fixed_macro_folder != "") {
    macro_folder = fixed_macro_folder;
    if (!endsWith(macro_folder, "/") && !endsWith(macro_folder, "\\")) macro_folder = macro_folder + "/";
    print("Merged CSV folder: using fixed_macro_folder (" + macro_folder + ")");
} else if (looks_valid) {
    print("Merged CSV folder: auto-detected from macro path (" + macro_folder + ")");
} else {
    macro_folder = output_folder;
    print("Merged CSV folder: could not auto-detect this macro's own path " +
          "(getInfo('macro.filepath') returned '" + macro_path + "'), so using " +
          "output_folder instead. Set fixed_macro_folder above to fix this.");
}

crops_folder = output_folder + "output_images/";
csv_folder   = output_folder + "analysis_csv/";
maps_folder  = output_folder + "cell_maps/";
masks_folder = output_folder + "cell_masks/";
skeleton_examples_folder = output_folder + "tagged_skeleton_examples/";
if (!File.exists(crops_folder)) File.makeDirectory(crops_folder);
if (!File.exists(csv_folder))   File.makeDirectory(csv_folder);
if (!File.exists(maps_folder))  File.makeDirectory(maps_folder);
if (!File.exists(masks_folder)) File.makeDirectory(masks_folder);
if (num_skeleton_examples > 0 && !File.exists(skeleton_examples_folder)) File.makeDirectory(skeleton_examples_folder);
reservoir_files = newArray(num_skeleton_examples);

// Single merged CSV across the whole run, saved directly next to this
// macro file, with the input folder's name (as an experiment tag) and
// the run's date/time baked into the filename. Columns match what
// data_helpers.R's load_morphology_data() produces when it merges the
// per-image CSVs post-hoc.
input_folder_trimmed = input_folder;
if (endsWith(input_folder_trimmed, "/") || endsWith(input_folder_trimmed, "\\"))
    input_folder_trimmed = substring(input_folder_trimmed, 0, lengthOf(input_folder_trimmed) - 1);
experiment_tag = replace(File.getName(input_folder_trimmed), "[^A-Za-z0-9]", "_");

getDateAndTime(run_yr, run_mo, run_dow, run_dom, run_hr, run_mn, run_sc, run_ms);
run_date_str = "" + run_yr + (run_mo + 1) + run_dom + run_hr + run_mn + run_sc;

merged_csv_path = macro_folder + "merged_morphology_data_" + experiment_tag + "_" + run_date_str + ".csv";
merged_header = "sample_id,image_id,cell_number,cell_name,area_um2,convexhull_area_um2,perimeter_um,convexhull_perimeter_um," +
                 "ramification_index,circularity,solidity,n_branches,n_junctions,tree_length_um," +
                 "avg_branch_length_um,max_branch_length_um,n_tips,n_triple_points,n_quadruple_points," +
                 "aspect_ratio,polarity_offset_um,centroid_x_um,centroid_y_um,source_file," +
                 "project,well_position,site,wavelength,image,ramification_index_2d";
File.saveString(merged_header + "\n", merged_csv_path);


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
print("Masks: " + masks_folder);
print("Merged: " + merged_csv_path);
if (num_skeleton_examples > 0) print("Skeleton examples: " + skeleton_examples_folder + " (" + minOf(skeleton_seen_count, num_skeleton_examples) + " of " + skeleton_seen_count + " seen)");
showMessage("Complete",
    "Processed " + ni + " image(s), " + grand_total_cells + " cells total.\n\n" +
    "Crops: output_images/ (<image_id>_cell<N>.tif)\n" +
    "Maps:  cell_maps/     (<image_id>_map.tif, numbered)\n" +
    "CSVs:  analysis_csv/ (R-pipeline ready)\n" +
    "Masks: cell_masks/    (<image_id>_cell<N>_mask.tif + mask_index_*.csv)\n" +
    "         -> feed cell_masks/ into the Channel Intensity Quantification\n" +
    "            pipeline to measure w1/w3/w4 within these same cell masks.\n\n" +
    "Merged CSV: saved next to this macro as\n" +
    "            merged_morphology_data_" + experiment_tag + "_" + run_date_str + ".csv\n" +
    "            (all images, all cells, same columns data_helpers.R adds\n" +
    "            during its own merge)\n\n" +
    "Skeleton examples: tagged_skeleton_examples/ (" + minOf(skeleton_seen_count, num_skeleton_examples) + " random of " + skeleton_seen_count + " total cells)");




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

    // --- merged-CSV metadata fields, mirroring data_helpers.R's
    //     parse_cell_name() / parse_filename_metadata() token rules ---
    sample_id = firstToken(base_name);
    image_id  = base_name;

    if (matches(base_name, ".*_[0-9]+$")) {
        us = lastIndexOf(base_name, "_");
        parsed_image = substring(base_name, us + 1, lengthOf(base_name));
        meta_rest = substring(base_name, 0, us);
    } else {
        parsed_image = "";
        meta_rest = base_name;
    }
    meta_toks = split(meta_rest, "_");

    idx = findFirstMatch(meta_toks, "[A-Z]{1,2}[0-9]+");
    if (idx >= 0) { parsed_well = meta_toks[idx]; meta_toks = removeAt(meta_toks, idx); } else parsed_well = "";

    idx = findFirstMatch(meta_toks, "(?i)w[0-9]+");
    if (idx >= 0) { parsed_wavelength = meta_toks[idx]; meta_toks = removeAt(meta_toks, idx); } else parsed_wavelength = "";

    idx = findFirstMatch(meta_toks, "(?i)s[0-9]+");
    if (idx >= 0) { parsed_site = meta_toks[idx]; meta_toks = removeAt(meta_toks, idx); } else parsed_site = "";

    parsed_project = joinUnderscore(meta_toks);

    header = "cell_name,area_um2,convexhull_area_um2,perimeter_um,convexhull_perimeter_um," +
             "ramification_index,circularity,solidity,n_branches,n_junctions,tree_length_um," +
             "avg_branch_length_um,max_branch_length_um,n_tips,n_triple_points,n_quadruple_points," +
             "aspect_ratio,polarity_offset_um,centroid_x_um,centroid_y_um";
    File.saveString(header + "\n", csv_path);

    mask_index_path = masks_folder + "mask_index_" + clean_base + "_" + date_str + time_str + ".csv";
    mask_index_header = "cell_number,cell_name,source_image,source_channel_token,mask_image," +
                         "crop_x,crop_y,crop_w,crop_h";
    File.saveString(mask_index_header + "\n", mask_index_path);

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

     
        setOption("BlackBackground", true);
        run("Create Mask");
        mask_name = base_name + "_cell" + (k + 1) + "_mask.tif";
        saveAs("TIFF", masks_folder + mask_name);
        close();
        selectWindow(crop_title);

        mask_row = d2s(k + 1, 0) + ",\"" + crop_name + "\",\"" + image_filename + "\",\"" +
                   channel_filter + "\",\"" + mask_name + "\"," +
                   d2s(cx, 0) + "," + d2s(cy, 0) + "," + d2s(cw, 0) + "," + d2s(ch, 0);
        File.append(mask_row, mask_index_path);
        // --- end v1.5 addition ---

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
        maybeSaveSkeletonExample(base_name + "_cell" + (k + 1), crop_title);
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

        ramification_2d = 0;
        if (area > 0) ramification_2d = perim / (2 * sqrt(PI * area));
        merged_row = "\"" + sample_id + "\",\"" + image_id + "\"," + (k + 1) + ",\"" + crop_name + "\"," +
              d2s(area,4) + "," + d2s(convex_area,4) + "," +
              d2s(perim,4) + "," + d2s(convex_perim,4) + "," +
              d2s(ramification,4) + "," + d2s(circ,4) + "," + d2s(solidity,4) + "," +
              n_branches + "," + n_junctions + "," + d2s(tree_len,4) + "," +
              d2s(avg_branch,4) + "," + d2s(max_branch,4) + "," +
              n_tips + "," + n_triple + "," + n_quad + "," +
              d2s(ar,4) + "," + d2s(polarity,4) + "," +
              d2s(centroid_x,4) + "," + d2s(centroid_y,4) + "," +
              "\"" + File.getName(csv_path) + "\"," +
              "\"" + parsed_project + "\",\"" + parsed_well + "\",\"" + parsed_site + "\",\"" + parsed_wavelength + "\",\"" + parsed_image + "\"," +
              d2s(ramification_2d,4);
        File.append(merged_row, merged_csv_path);

        
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


function maybeSaveSkeletonExample(label, base_window_title) {
    // Reservoir sampling (Algorithm R): keeps a uniformly random sample of
    // size num_skeleton_examples from the whole stream of cells processed,
    // in a single pass, without knowing the total count ahead of time.
    if (num_skeleton_examples <= 0) return;
    if (!isOpen("Tagged skeleton")) return;

    skeleton_seen_count++;
    fname = "example_tagged_skeleton_" + label + ".tif";

    if (skeleton_seen_count <= num_skeleton_examples) {
        // reservoir not yet full: always keep
        slot = skeleton_seen_count - 1;
        saveSkeletonOverlay(base_window_title, skeleton_examples_folder + fname);
        reservoir_files[slot] = fname;
    } else {
        // reservoir full: replace a random existing slot with
        // probability num_skeleton_examples / skeleton_seen_count
        r = floor(random() * skeleton_seen_count);
        if (r < num_skeleton_examples) {
            old_fname = reservoir_files[r];
            if (File.exists(skeleton_examples_folder + old_fname)) File.delete(skeleton_examples_folder + old_fname);
            saveSkeletonOverlay(base_window_title, skeleton_examples_folder + fname);
            reservoir_files[r] = fname;
        }
    }
}

function saveSkeletonOverlay(base_window_title, save_path) {
    // Composites the color-coded "Tagged skeleton" on top of this cell's
    // cropped image, so the example shows the traced skeleton directly on
    // the real cell instead of alone on a black background.
    // Both inputs are duplicated first so the originals (still needed
    // elsewhere / closed by the caller) are left untouched.
    selectWindow(base_window_title);
    run("Select None");
    run("Duplicate...", "title=overlay_base_tmp");
    if (bitDepth() != 24) run("RGB Color");

    selectWindow("Tagged skeleton");
    run("Select None");
    run("Duplicate...", "title=overlay_skel_tmp");
    if (bitDepth() != 24) run("RGB Color");

    // Brightest-pixel-wins blend: colored skeleton pixels are brighter than
    // the black background they sit on, so they show through on top of the
    // cell wherever they occur, while true background leaves the cell as-is.
    imageCalculator("Max create", "overlay_base_tmp", "overlay_skel_tmp");
    overlay_title = getTitle();

    saveAs("TIFF", save_path);
    close(overlay_title);
    close("overlay_base_tmp");
    close("overlay_skel_tmp");
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


function firstToken(s) {
    // First underscore-delimited token, mirroring data_helpers.R's
    // sample_id <- str_extract(cell_name, "^[^_]+")
    idx = indexOf(s, "_");
    if (idx < 0) return s;
    return substring(s, 0, idx);
}

function findFirstMatch(toks, patternRegex) {
    // Index of the first array element that fully matches patternRegex,
    // or -1 if none do. Mirrors data_helpers.R's take()/str_which() step.
    for (i = 0; i < toks.length; i++) {
        if (matches(toks[i], patternRegex)) return i;
    }
    return -1;
}

function removeAt(toks, idx) {
    // New array with the element at idx dropped.
    out = newArray(toks.length - 1);
    j = 0;
    for (i = 0; i < toks.length; i++) {
        if (i != idx) { out[j] = toks[i]; j++; }
    }
    return out;
}

function joinUnderscore(toks) {
    // Rejoins whatever tokens are left (after well/wavelength/site were
    // pulled out) into the "project" field, same as data_helpers.R's
    // paste(toks, collapse = "_").
    s = "";
    for (i = 0; i < toks.length; i++) {
        if (i > 0) s = s + "_";
        s = s + toks[i];
    }
    return s;
}
