// v1.5: exports the per cell mask plus a mask_index CSV linking each mask
//       to its source image and crop rectangle, for the Channel Intensity
//       Quantification pipeline.

// v1.6: saves a random sample of tagged skeleton images instead of just
//       the first one, into tagged_skeleton_examples/.

// v1.7: skeleton examples composite onto the cropped cell image. Output
//       folder auto detects from the macro file location. Per image CSVs
//       also merge into one merged_morphology_data.csv with the same
//       derived columns data_helpers.R adds. fixed_output_folder and
//       fixed_macro_folder let you override auto detection when needed.

// v1.8: A whole image outlier check against a running batch baseline, an
//       absolute check for essentially empty images, an optional per cell
//       shape check (off by default), brightness/contrast inspection, an
//       auto brightness fix for zero candidate images, and automatic
//       removal of tiny noise particles. Any flagged image can also turn
//       off further flagging for the rest of the run. Review dialogs are
//       positioned so they do not cover the image, and the Results table
//       auto closes before they appear.
//       Key tunables: outlier_min_baseline, outlier_z_threshold,
//       outlier_mad_floor_fraction, empty_image_max_count,
//       empty_image_max_area_fraction, refined_area_fraction,
//       refined_min_ramification, auto_exclude_area_fraction,
//       auto_brightness_saturated_pct, enable_per_cell_shape_check.

// v1.9: Merged CSV now saves into the auto-generated analysis_csv/ folder
//       (next to the per-image CSVs) instead of the macro file's own
//       folder. fixed_macro_folder is removed


var threshold_algorithm, min_cell_size, max_cell_size;
var subtract_bg_radius, smooth_radius, close_gaps, fill_detection_holes;
var split_touching, split_dynamic, mask_grow, padding, fallback_um_per_px;
var input_folder, output_folder, crops_folder, csv_folder, maps_folder, masks_folder;
var images, ni, grand_total_cells, preset, accept_all;
var accept_all_batch, border_margin;
var num_skeleton_examples, skeleton_examples_folder, skeleton_seen_count, reservoir_files;
var merged_csv_path;
var baseline_counts, baseline_areas, outlier_min_baseline, outlier_z_threshold;
var outlier_mad_floor_fraction, outlier_checks_disabled;
var refined_area_fraction, refined_min_ramification, enable_per_cell_shape_check;
var auto_exclude_area_fraction;
var auto_brightness_fix_enabled, auto_brightness_saturated_pct;
var empty_image_max_count, empty_image_max_area_fraction;

outlier_min_baseline = 8;
outlier_z_threshold  = 5.0;
outlier_mad_floor_fraction = 0.25;
outlier_checks_disabled = false;
baseline_counts = newArray();
baseline_areas  = newArray();
auto_brightness_fix_enabled = false;
auto_brightness_saturated_pct = 0.35;
empty_image_max_count         = 3;
empty_image_max_area_fraction = 0.5;
auto_exclude_area_fraction = 0.15;
enable_per_cell_shape_check = false;
refined_area_fraction    = 0.5;
refined_min_ramification = 1.15;

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


fixed_output_folder = "";

applyPreset(preset);      


if (preset == "Small cells")            fallback_um_per_px = 4 * (184.65 / 512);
else if (preset == "ImageXpress") fallback_um_per_px = 0.65;  
else                                        fallback_um_per_px = 184.65 / 512;


input_folder  = getDirectory("Choose the input folder (original images)");

macro_path = getInfo("macro.filepath");
looks_valid = (macro_path != "") && endsWith(toLowerCase(macro_path), ".ijm");
if (looks_valid) macro_folder = File.getDirectory(macro_path);

if (fixed_output_folder != "") {
    output_folder = fixed_output_folder;
    if (!endsWith(output_folder, "/") && !endsWith(output_folder, "\\")) output_folder = output_folder + "/";
} else if (looks_valid) {
    output_folder = macro_folder;
} else {
    output_folder = getDirectory("Choose the output folder (results go here)");
}
if (!File.exists(output_folder)) File.makeDirectory(output_folder);

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

input_folder_trimmed = input_folder;
if (endsWith(input_folder_trimmed, "/") || endsWith(input_folder_trimmed, "\\"))
    input_folder_trimmed = substring(input_folder_trimmed, 0, lengthOf(input_folder_trimmed) - 1);
experiment_tag = replace(File.getName(input_folder_trimmed), "[^A-Za-z0-9]", "_");

getDateAndTime(run_yr, run_mo, run_dow, run_dom, run_hr, run_mn, run_sc, run_ms);
run_date_str = "" + run_yr + (run_mo + 1) + run_dom + run_hr + run_mn + run_sc;

merged_csv_path = csv_folder + "merged_morphology_data_" + experiment_tag + "_" + run_date_str + ".csv";
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
    "Crops: output_images/\n" +
    "Maps:  cell_maps/\n" +
    "CSVs:  analysis_csv/\n" +
    "Masks: cell_masks/ - use for channel intensity quantification\n" +
    "Merged CSV: saved in parent directory\n" +
    "Skeleton examples: tagged_skeleton_examples/ (" + minOf(skeleton_seen_count, num_skeleton_examples) + " of " + skeleton_seen_count + ")");




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
    n_exported = 0;

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

            layoutForReview(original_title);
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

        pre_count = roiManager("count");
        run("Set Measurements...", "area redirect=None decimal=4");
        run("Clear Results");
        run("Analyze Particles...", "size=0-Infinity pixel show=Nothing display add");
        resetThreshold();
        n_particles = roiManager("count") - pre_count;

        if (n_particles == 0) {
            area = 0;
            refined_idx = -1;
        } else {
            largest_area = -1;
            largest_p = 0;
            for (p = 0; p < n_particles; p++) {
                a = getResult("Area", p);
                if (a > largest_area) { largest_area = a; largest_p = p; }
            }
            area = largest_area;

            roiManager("select", pre_count + largest_p);
            roiManager("add");
            kept_idx = roiManager("count") - 1;
            to_delete = newArray(n_particles);
            for (p = 0; p < n_particles; p++) to_delete[p] = pre_count + p;
            roiManager("select", to_delete);
            roiManager("delete");
            refined_idx = kept_idx - n_particles;
        }

        if (refined_idx == -1 || area < auto_exclude_area_fraction * min_cell_size) {
            if (refined_idx != -1) { roiManager("select", refined_idx); roiManager("delete"); }
            close(crop_title);
            continue;
        }

        roiManager("select", refined_idx);
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

        roiManager("select", refined_idx);
        run("Convex Hull");
        run("Measure");
        convex_area  = getResult("Area", 1);
        convex_perim = getResult("Perim.", 1);
        if (convex_perim > 0) ramification = perim / convex_perim; else ramification = 0;

        // --- pause on cells that look too small/too convex to be a
        //     confident ramified-microglia call, when accept_all(_batch)
        //     would otherwise wave them through with no human ever seeing
        //     them. A "yes, keep it" here is just as valid an answer as
        //     "no" - this is a prompt, not a filter. OFF by default (see
        //     enable_per_cell_shape_check near the top) whole-image
        //     flagging below still runs regardless.
        is_suspicious = enable_per_cell_shape_check &&
                        ((area < refined_area_fraction * min_cell_size) ||
                         (ramification < refined_min_ramification));
        keep_cell = true;
        if (is_suspicious && (accept_all || accept_all_batch)) {
            keep_cell = confirmSuspiciousCell(crop_title, refined_idx, image_filename,
                                               k + 1, area, ramification);
        }
        if (!keep_cell) {
            roiManager("select", refined_idx);
            roiManager("delete");
            close(crop_title);
            continue;
        }
        n_exported++;

        // export the mask from the SAVED refined selection (not
        //     whatever selection happens to be active after the Convex
        //     Hull step above, that would be the hull, not the true shape)
        newImage("mask_export_tmp", "8-bit black", cw, ch, 1);
        roiManager("select", refined_idx);
        setColor(255);
        fill();
        run("Select None");
        mask_name = base_name + "_cell" + (k + 1) + "_mask.tif";
        saveAs("TIFF", masks_folder + mask_name);
        close();
        roiManager("select", refined_idx);
        roiManager("delete");
        selectWindow(crop_title);

        mask_row = d2s(k + 1, 0) + ",\"" + crop_name + "\",\"" + image_filename + "\",\"" +
                   channel_filter + "\",\"" + mask_name + "\"," +
                   d2s(cx, 0) + "," + d2s(cy, 0) + "," + d2s(cw, 0) + "," + d2s(ch, 0);
        File.append(mask_row, mask_index_path);
     

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
    grand_total_cells += n_exported;
    if (n_exported != napproved) {
        print("  " + (napproved - n_exported) + " approved cell(s) dropped as noise/debris after re-thresholding, or skipped via prompt");
    }
    print("  CSV: " + csv_path);
}


function maybeSaveSkeletonExample(label, base_window_title) {

    if (num_skeleton_examples <= 0) return;
    if (!isOpen("Tagged skeleton")) return;

    skeleton_seen_count++;
    fname = "example_tagged_skeleton_" + label + ".tif";

    if (skeleton_seen_count <= num_skeleton_examples) {
    
        slot = skeleton_seen_count - 1;
        saveSkeletonOverlay(base_window_title, skeleton_examples_folder + fname);
        reservoir_files[slot] = fname;
    } else {
      
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
    // Both inputs are duplicated first so the originals
    selectWindow(base_window_title);
    run("Select None");
    run("Duplicate...", "title=overlay_base_tmp");
    if (bitDepth() != 24) run("RGB Color");

    selectWindow("Tagged skeleton");
    run("Select None");
    run("Duplicate...", "title=overlay_skel_tmp");
    if (bitDepth() != 24) run("RGB Color");

  
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


function inspectBrightnessContrast(image_title) {
    // Adjusts, then unconditionally bakes whatever range you land on into
    // the actual pixel values (not just the display) once you click OK
    // no need to separately click Apply on the B&C panel itself. Returns
    // [min, max] that was applied, so a caller can save/replay it exactly.
    selectWindow(image_title);
    run("Brightness/Contrast...");
    if (isOpen("B&C")) { selectWindow("B&C"); setLocation(20, 460); }
    selectWindow(image_title);
    waitForUser("Inspect brightness/contrast",
        "Manually adjust brightness and contrast, then press OK here.");
    if (isOpen("B&C")) { selectWindow("B&C"); run("Close"); }
    selectWindow(image_title);
    getMinAndMax(applied_min, applied_max);
    setMinAndMax(applied_min, applied_max);
    run("Apply LUT");
    result = newArray(2);
    result[0] = applied_min;
    result[1] = applied_max;
    return result;
}

function layoutForReview(image_title) {
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    if (isOpen(image_title)) {
        selectWindow(image_title);
        setLocation(20, 60);
    }
    dialog_x = screenWidth - 460;
    if (dialog_x < 440) dialog_x = 440;
    Dialog.setLocation(dialog_x, 60);
}

function tuneDetection(original_title, width, height, image_filename) {
    if (accept_all_batch) {
        detectCells(original_title, width, height);
        stats = candidateStats();
        n_candidates = stats[0];
        med_area = stats[1];

        // If auto-brightness-fix was enabled (see below) and THIS image
        // also found zero candidates
        if (n_candidates == 0 && auto_brightness_fix_enabled) {
            selectWindow(original_title);
            run("Enhance Contrast...", "saturated=" + auto_brightness_saturated_pct + " normalize");
            detectCells(original_title, width, height);
            stats = candidateStats();
            n_candidates = stats[0];
            med_area = stats[1];
        }

        reason = checkImageOutlier(n_candidates, med_area);

        if (reason == "") {
            recordBaseline(n_candidates, med_area);
            return true;
        }

        // flagged: pause just for this one image, batch mode resumes after
        selectWindow(original_title);
        Overlay.remove;
        for (r = 0; r < n_candidates; r++) { roiManager("select", r); Overlay.addSelection("magenta"); }
        Overlay.show;
        run("Select None");

        layoutForReview(original_title);
        Dialog.create("Unusual image - " + image_filename);
        Dialog.addMessage(image_filename);
        Dialog.addMessage(reason);
        if (n_candidates == 0) {
            Dialog.addChoice("What do you want to do with this image?",
                newArray("Skip it entirely (treat as 0 cells)", "Inspect brightness/contrast, then decide",
                         "Adjust settings and re-detect", "Stop flagging for the rest of this run"),
                "Inspect brightness/contrast, then decide");
        } else {
            Dialog.addMessage("Magenta outlines = " + n_candidates + " candidate(s) found. " +
                               "'Accept ALL for ALL remaining\nimages' stays on for every other image.");
            Dialog.addChoice("What do you want to do with this image?",
                newArray("Review it manually (KEEP/SKIP/SPLIT)", "Accept it anyway (keep all candidates)",
                         "Skip it entirely (treat as 0 cells)", "Inspect brightness/contrast, then decide",
                         "Adjust settings and re-detect", "Stop flagging for the rest of this run"),
                "Review it manually (KEEP/SKIP/SPLIT)");
        }
        Dialog.show();
        action = Dialog.getChoice();
        Overlay.remove;

        if (action == "Stop flagging for the rest of this run") {
            outlier_checks_disabled = true;
            print("  Outlier flagging disabled for the rest of this run - Accept ALL will no longer interrupt.");
            recordBaseline(n_candidates, med_area);
            return true;
        }

        
        if (action == "Accept it anyway (keep all candidates)") {
            recordBaseline(n_candidates, med_area);
            return true;
        }

        if (action == "Inspect brightness/contrast, then decide") {
            was_zero = (n_candidates == 0);
            bc_range = inspectBrightnessContrast(original_title);
            if (was_zero) {
                layoutForReview(original_title);
                Dialog.create("Apply this fix automatically?");
                Dialog.addMessage("Apply an automatic brightness fix to all images that also find zero candidates?");
                Dialog.addCheckbox("Yes, auto-fix future zero-candidate images", true);
                Dialog.show();
                if (Dialog.getCheckbox()) {
                    auto_brightness_fix_enabled = true;
                    print("  Auto-brightness-fix enabled for future zero-candidate images.");
                }
            }
            return tuneDetection(original_title, width, height, image_filename);
        }
        if (action == "Adjust settings and re-detect") {
            advancedDialog();
            return tuneDetection(original_title, width, height, image_filename);
        }
        if (action == "Skip it entirely (treat as 0 cells)") {
            roiManager("reset");
            return true;
        }
        if (action == "Review it manually (KEEP/SKIP/SPLIT)") {
            return false;
        }
        return true;   // "Accept it anyway" fallback 
    }

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

        layoutForReview(original_title);
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
    // This image was looked at directly (by the person reading this right now), so its
    // stats are as trustworthy a baseline sample as any feed them in
    // before "Accept ALL for ALL remaining images" starts relying on it.
    stats = candidateStats();
    recordBaseline(stats[0], stats[1]);

    Overlay.remove;
    return accept_all;
}

function confirmSuspiciousCell(crop_title, refined_idx, image_filename, cell_num, area, ramification) {
   
    selectWindow(crop_title);
    Overlay.remove;
    roiManager("select", refined_idx);
    Overlay.addSelection("magenta");
    Overlay.show;
    run("Select None");

    layoutForReview(crop_title);
    Dialog.create("Unusual cell shape - cell " + cell_num);
    Dialog.addMessage(image_filename + "  (cell " + cell_num + ")");
    Dialog.addMessage("Area " + d2s(area, 0) + " px, ramification " + d2s(ramification, 2) +
                       " (1.0 = perfectly round/convex, higher = more branched).\n\n" +
                       "This is smaller/rounder than typical for a ramified microglia - could be\n" +
                       "debris or an imaging artifact, OR a genuine ameboid/activated microglia.\n" +
                       "Check the outlined shape (magenta) and decide.");
    Dialog.addChoice("Keep this as a cell?",
        newArray("SKIP - not a real cell", "KEEP - this is a real cell", "Let me look closer first (brightness/contrast)"),
        "SKIP - not a real cell");
    Dialog.show();
    choice = Dialog.getChoice();
    Overlay.remove;
    if (startsWith(choice, "Let me look closer")) {
        inspectBrightnessContrast(crop_title);
        return confirmSuspiciousCell(crop_title, refined_idx, image_filename, cell_num, area, ramification);
    }
    return startsWith(choice, "KEEP");
}

function candidateStats() {
    n = roiManager("count");
    result = newArray(2);
    result[0] = n;
    if (n == 0) { result[1] = 0; return result; }
    areas = newArray(n);
    for (r = 0; r < n; r++) {
        roiManager("select", r);
        getStatistics(a);
        areas[r] = a;
    }
    roiManager("deselect");
    result[1] = computeMedian(areas);
    return result;
}

function computeMedian(arr) {
    n = arr.length;
    if (n == 0) return 0;
    tmp = Array.copy(arr);
    Array.sort(tmp);
    mid = floor(n / 2);
    if (n % 2 == 0) return (tmp[mid - 1] + tmp[mid]) / 2;
    return tmp[mid];
}

function computeMAD(arr, med) {
    n = arr.length;
    if (n == 0) return 0;
    devs = newArray(n);
    for (i = 0; i < n; i++) devs[i] = abs(arr[i] - med);
    return computeMedian(devs);
}

function recordBaseline(n_candidates, med_area) {
    baseline_counts[baseline_counts.length] = n_candidates;
    if (n_candidates > 0) baseline_areas[baseline_areas.length] = med_area;
}

function checkImageOutlier(n_candidates, med_area) {
    if (outlier_checks_disabled) return "";
    reason = "";

    // --- absolute check: does this image look essentially empty? ---
    // Unlike the relative checks below, this doesn't need any baseline -
    // it fires from image 1, and can't be diluted by other sparse images
    // getting folded into the baseline first.
    if (n_candidates <= empty_image_max_count) {
        if (n_candidates == 0) {
            reason = reason + "No candidates found at all in this image.\n";
        } else if (med_area < empty_image_max_area_fraction * min_cell_size) {
            reason = reason + "Only " + n_candidates + " candidate(s) found, and they're small " +
                     "(median " + d2s(med_area, 0) + " px, vs a typical minimum around " +
                     d2s(min_cell_size, 0) + " px) - this image may be essentially empty " +
                     "(wrong focal plane, no real cells, or pure noise/debris).\n";
        }
    }

    if (baseline_counts.length >= outlier_min_baseline) {
        med_count = computeMedian(baseline_counts);
        mad_count = computeMAD(baseline_counts, med_count);
        mad_count_floor = med_count * outlier_mad_floor_fraction;
        if (mad_count < mad_count_floor) mad_count = mad_count_floor;
        if (mad_count > 0) {
            z_count = 0.6745 * (n_candidates - med_count) / mad_count;
            if (abs(z_count) > outlier_z_threshold) {
                reason = reason + "Cell count " + n_candidates + " is unusual (batch median so far: " +
                         d2s(med_count, 1) + ").\n";
            }
        }
    }

    if (baseline_areas.length >= outlier_min_baseline && n_candidates > 0) {
        med_area_base = computeMedian(baseline_areas);
        mad_area_base = computeMAD(baseline_areas, med_area_base);
        mad_area_floor = med_area_base * outlier_mad_floor_fraction;
        if (mad_area_base < mad_area_floor) mad_area_base = mad_area_floor;
        if (mad_area_base > 0) {
            z_area = 0.6745 * (med_area - med_area_base) / mad_area_base;
            if (z_area < -outlier_z_threshold) {
                reason = reason + "Median candidate size " + d2s(med_area, 0) + " px is much smaller " +
                         "than the batch median so far (" + d2s(med_area_base, 0) + " px) - possible " +
                         "noise/debris instead of real cells (e.g. a wrong-focal-plane image).\n";
            }
        }
    }

    return reason;
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
    s = "";
    for (i = 0; i < toks.length; i++) {
        if (i > 0) s = s + "_";
        s = s + toks[i];
    }
    return s;
}
