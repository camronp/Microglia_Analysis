// v1.2 Channel Intensity Quantification
//
// Companion pipeline to v1_5_ImageXpress_Microglia_Segment_Analysis.ijm.

//

//
// v1.1 added:
//   - provides more example images of the mask it uses and the overlay onto fluorescence


var masks_folder, channel_input_folder, output_folder, results_csv_path, overlay_folder;
var channel_list, missing_log, overlay_counts, overlay_examples_per_channel;
var do_lcd_classification, lcd_channel_token, lcd_stat_key, lcd_stat_label, lcd_threshold, lcd_direction;

stat_labels = newArray("Mean intensity", "Median intensity", "Min intensity", "Max intensity",
                       "Std. dev. intensity", "Integrated density", "Raw integrated density", "Area (px)");
stat_keys   = newArray("mean", "median", "min", "max", "stddev", "integrated", "rawintegrated", "area");
direction_choices = newArray("Below threshold = unseeded", "Below threshold = seeded");

Dialog.create("Channel Intensity Quantification");
Dialog.addMessage("Applies the cell masks saved by the microglia segmentation\n" +
                   "pipeline to other channels, and measures intensity inside\n" +
                   "each cell's mask.");
Dialog.addString("Channels to quantify (comma-separated tokens, e.g. w1,w3,w4)", "w1,w3,w4");
Dialog.addCheckbox("Save example mask-overlay images (visual QC)", true);
Dialog.addNumber("Number of example overlays to save per channel", 10);
Dialog.addMessage("--- Optional: classify cells as seeded / unseeded ---");
Dialog.addCheckbox("Classify cells using one channel (e.g. LCD Fibrils)?", false);
Dialog.addString("That channel's token (e.g. w3)", "");
Dialog.addChoice("Statistic to threshold on", stat_labels, stat_labels[0]);
Dialog.addNumber("Threshold value", 200);
Dialog.addChoice("Classification direction", direction_choices, direction_choices[0]);
Dialog.show();

channel_field                = Dialog.getString();
save_overlay_examples        = Dialog.getCheckbox();
overlay_examples_per_channel = round(Dialog.getNumber());
do_lcd_classification         = Dialog.getCheckbox();
lcd_channel_token             = String.trim(Dialog.getString());
lcd_stat_label                = Dialog.getChoice();
lcd_threshold                 = Dialog.getNumber();
lcd_direction                 = Dialog.getChoice();

if (do_lcd_classification && lcd_channel_token == "") {
    showMessage("Missing channel token",
        "You checked 'Classify cells...' but didn't enter a channel token - " +
        "classification will be skipped.");
    do_lcd_classification = false;
}
lcd_stat_key = stat_keys[indexOfString(stat_labels, lcd_stat_label)];

masks_folder = getDirectory("Choose the cell_masks folder (output of the segmentation pipeline)");
channel_input_folder = getDirectory("Choose the folder containing the OTHER channel images (w1/w3/w4/...)");
output_folder = getDirectory("Choose the OUTPUT folder for the intensity results CSV");

// --- parse the requested channel token list ---
channel_list = split(channel_field, ",");
for (i = 0; i < channel_list.length; i++) channel_list[i] = String.trim(channel_list[i]);

overlay_folder = output_folder + "mask_overlay_examples/";
if (save_overlay_examples && !File.exists(overlay_folder)) File.makeDirectory(overlay_folder);
overlay_counts = newArray(channel_list.length);   // one running count per channel_list entry, starts at 0

print("\\Clear");
print("=== Channel Intensity Quantification ===");
print("Masks folder  : " + masks_folder);
print("Channel folder: " + channel_input_folder);
print("Output folder : " + output_folder);
print("Channels      : " + channel_field);
if (do_lcd_classification) {
    print("Seeded/unseeded: channel '" + lcd_channel_token + "', " + lcd_stat_label +
          " " + lcd_direction + " " + d2s(lcd_threshold, 4));
}
print("");

// --- find all mask_index CSVs (one per source image processed upstream) ---
mask_files = getFileList(masks_folder);
index_files = newArray();
nidx = 0;
for (f = 0; f < mask_files.length; f++) {
    if (startsWith(mask_files[f], "mask_index_") && endsWith(mask_files[f], ".csv")) {
        index_files[nidx] = mask_files[f];
        nidx++;
    }
}

if (nidx == 0) {
    showMessage("No mask index found",
        "No 'mask_index_*.csv' files were found in:\n" + masks_folder +
        "\n\nMake sure you selected the cell_masks/ folder produced by the\n" +
        "v1.5 segmentation pipeline (older v1.4 outputs don't have these).");
    exit;
}
print("Found " + nidx + " mask_index file(s)");

// --- prepare the combined output CSV (long/tidy format: 1 row per cell x channel) ---
getDateAndTime(yr, mo, dow, dom, hr, mn, sc, ms);
ts = IJ.pad(yr,4) + IJ.pad(mo+1,2) + IJ.pad(dom,2) + "_" +
     IJ.pad(hr,2) + IJ.pad(mn,2) + IJ.pad(sc,2);
results_csv_path = output_folder + "channel_intensity_results_" + ts + ".csv";
results_header = "cell_name,source_image,channel,mean_intensity,min_intensity,max_intensity," +
                  "median_intensity,stddev_intensity,integrated_density,raw_integrated_density,area_px," +
                  "lcd_channel,lcd_statistic,lcd_value,seeded_status";
File.saveString(results_header + "\n", results_csv_path);

missing_log = newArray();
nmissing = 0;
n_measurements = 0;
n_seeded = 0;
n_unseeded = 0;

setBatchMode(true);
total_rows_seen = 0;
total_rows_malformed = 0;

for (f = 0; f < nidx; f++) {
    content = File.openAsString(masks_folder + index_files[f]);
    lines = split(content, "\n");
    rows_this_file = 0;
    for (li = 1; li < lines.length; li++) {   // li=0 is the header row, skip it
        line = String.trim(lines[li]);
        if (lengthOf(line) == 0) continue;
        rows_this_file++;
        total_rows_seen++;
        fields = parseCSVLine(line);
        if (fields.length < 9) {
            total_rows_malformed++;
            snippet = line;
            if (lengthOf(snippet) > 80) snippet = substring(snippet, 0, 80) + "...";
            logMissing("(unknown)", index_files[f],
                "(all channels)",
                "malformed row (" + fields.length + " field(s), expected 9): " + snippet);
            continue;
        }

        cell_number         = fields[0];
        cell_name           = fields[1];
        source_image        = fields[2];
        source_channel_token = fields[3];
        mask_image          = fields[4];
        crop_x = parseFloat(fields[5]);
        crop_y = parseFloat(fields[6]);
        crop_w = parseFloat(fields[7]);
        crop_h = parseFloat(fields[8]);

        // --- optional: classify this cell as seeded/unseeded from one channel ---
        cell_lcd_value = NaN;
        cell_seeded_status = "";
        if (do_lcd_classification) {
            lcd_derived_name = replace(source_image, source_channel_token, lcd_channel_token);
            if (lcd_derived_name == source_image) {
                logMissing(cell_name, source_image, "(LCD:" + lcd_channel_token + ")",
                           "LCD channel token not found in filename");
            } else if (!File.exists(channel_input_folder + lcd_derived_name)) {
                logMissing(cell_name, source_image, "(LCD:" + lcd_channel_token + ")",
                           "LCD file not found: " + lcd_derived_name);
            } else if (!File.exists(masks_folder + mask_image)) {
                logMissing(cell_name, source_image, "(LCD:" + lcd_channel_token + ")",
                           "mask file missing: " + mask_image);
            } else {
                open(channel_input_folder + lcd_derived_name);
                lcd_title = getTitle();
                makeRectangle(crop_x, crop_y, crop_w, crop_h);
                run("Crop");
                run("Select None");

                open(masks_folder + mask_image);
                lcd_mask_title = getTitle();
                setThreshold(128, 255);
                run("Create Selection");
                resetThreshold();
                close(lcd_mask_title);

                selectWindow(lcd_title);
                run("Restore Selection");
                if (selectionType() == -1) {
                    logMissing(cell_name, source_image, "(LCD:" + lcd_channel_token + ")",
                               "empty mask - no cell pixels found");
                } else {
                    run("Clear Results");
                    run("Set Measurements...",
                        "area mean standard min median integrated redirect=None decimal=4");
                    run("Measure");
                    cell_lcd_value = pick_stat(lcd_stat_key,
                        getResult("Mean", 0), getResult("Min", 0), getResult("Max", 0),
                        getResult("Median", 0), getResult("StdDev", 0),
                        getResult("IntDen", 0), getResult("RawIntDen", 0), getResult("Area", 0));

                    below = cell_lcd_value < lcd_threshold;
                    if (lcd_direction == "Below threshold = unseeded") {
                        if (below) cell_seeded_status = "unseeded"; else cell_seeded_status = "seeded";
                    } else {
                        if (below) cell_seeded_status = "seeded"; else cell_seeded_status = "unseeded";
                    }
                    if (cell_seeded_status == "seeded") n_seeded++; else n_unseeded++;
                }
                close(lcd_title);
            }
        }

        for (c = 0; c < channel_list.length; c++) {
            target_token = channel_list[c];
            if (target_token == "") continue;

            derived_name = replace(source_image, source_channel_token, target_token);
            if (derived_name == source_image) {
                logMissing(cell_name, source_image, target_token, "channel token not found in filename");
                continue;
            }
            if (!File.exists(channel_input_folder + derived_name)) {
                logMissing(cell_name, source_image, target_token, "file not found: " + derived_name);
                continue;
            }
            if (!File.exists(masks_folder + mask_image)) {
                logMissing(cell_name, source_image, target_token, "mask file missing: " + mask_image);
                continue;
            }

            // --- crop the target channel image to the same rectangle used for the mask ---
            open(channel_input_folder + derived_name);
            chan_title = getTitle();
            makeRectangle(crop_x, crop_y, crop_w, crop_h);
            run("Crop");
            run("Select None");

            // --- load the mask, threshold it, transfer the selection onto the crop ---
            open(masks_folder + mask_image);
            mask_title = getTitle();
            setThreshold(128, 255);
            run("Create Selection");
            resetThreshold();
            close(mask_title);

            selectWindow(chan_title);
            run("Restore Selection");
            if (selectionType() == -1) {
                logMissing(cell_name, source_image, target_token, "empty mask - no cell pixels found");
                close(chan_title);
                continue;
            }

            run("Clear Results");
            run("Set Measurements...",
                "area mean standard min median integrated redirect=None decimal=4");
            run("Measure");
            mean_val   = getResult("Mean", 0);
            min_val    = getResult("Min", 0);
            max_val    = getResult("Max", 0);
            median_val = getResult("Median", 0);
            stddev_val = getResult("StdDev", 0);
            intden_val = getResult("IntDen", 0);
            rawintden  = getResult("RawIntDen", 0);
            area_px    = getResult("Area", 0);

            row = "\"" + cell_name + "\",\"" + source_image + "\",\"" + target_token + "\"," +
                  d2s(mean_val,4) + "," + d2s(min_val,4) + "," + d2s(max_val,4) + "," +
                  d2s(median_val,4) + "," + d2s(stddev_val,4) + "," +
                  d2s(intden_val,4) + "," + d2s(rawintden,4) + "," + d2s(area_px,4) + "," +
                  lcd_row_channel(do_lcd_classification, lcd_channel_token) + "," +
                  lcd_row_stat(do_lcd_classification, lcd_stat_label) + "," +
                  lcd_row_value(cell_lcd_value) + "," +
                  "\"" + cell_seeded_status + "\"";
            File.append(row, results_csv_path);
            n_measurements++;

            // --- save example overlays per channel, for a visual sanity check ---
            if (save_overlay_examples && overlay_counts[c] < overlay_examples_per_channel) {
                run("Select None");
                run("Duplicate...", "title=overlay_base");
                run("Restore Selection");
                run("Enhance Contrast...", "saturated=0.35 normalize");
                run("RGB Color");
                setForegroundColor(255, 255, 0);
                setLineWidth(2);
                run("Draw");
                setColor(255, 255, 0);
                setFont("SansSerif", 12, "bold");
                drawString(cell_name + " (" + target_token + ")", 4, 16);
                overlay_name = "example_overlay_" + target_token + "_" +
                               IJ.pad(overlay_counts[c] + 1, 2) + ".png";
                saveAs("PNG", overlay_folder + overlay_name);
                close();
                selectWindow(chan_title);
                overlay_counts[c] = overlay_counts[c] + 1;
            }

            close(chan_title);
        }
    }
    print(index_files[f] + ": " + rows_this_file + " cell row(s)");
}
setBatchMode(false);

print("");
print("=== DONE ===");
print("Total data rows read : " + total_rows_seen);
print("Malformed rows       : " + total_rows_malformed);
print("Measurements written: " + n_measurements);
print("Results CSV: " + results_csv_path);
if (do_lcd_classification) {
    print("Seeded: " + n_seeded + "   Unseeded: " + n_unseeded);
}
if (save_overlay_examples) {
    print("Mask overlay examples: " + overlay_folder);
}

if (missing_log.length > 0) {
    print("");
    print("--- " + missing_log.length + " skipped (see details) ---");
    for (i = 0; i < missing_log.length; i++) print("  " + missing_log[i]);
}

summary_msg = "Wrote " + n_measurements + " measurement(s) to:\n" + results_csv_path;
if (do_lcd_classification) {
    summary_msg = summary_msg + "\n\nSeeded: " + n_seeded + "   Unseeded: " + n_unseeded;
}
if (save_overlay_examples) {
    summary_msg = summary_msg + "\n\nExample mask overlays saved to:\n" + overlay_folder;
}
if (missing_log.length > 0) {
    summary_msg = summary_msg + "\n\n" + missing_log.length +
        " cell/channel combination(s) were skipped - see the Log window for details.";
}
showMessage("Complete", summary_msg);


function logMissing(cell_name, source_image, target_token, reason) {
    missing_log[missing_log.length] = cell_name + " (" + source_image + " -> " + target_token + "): " + reason;
}

// These three build CSV fields as proper strings from the start (never a
// raw number leading a "+" chain) - see the note in v1.5's mask_row fix:
// ImageJ macro can silently turn "number + string + string..." into NaN
// if the chain doesn't start out as a string.
function lcd_row_channel(active, token) {
    if (!active) return "\"\"";
    return "\"" + token + "\"";
}
function lcd_row_stat(active, label) {
    if (!active) return "\"\"";
    return "\"" + label + "\"";
}
function lcd_row_value(value) {
    if (isNaN(value)) return "";
    return d2s(value, 4);
}

// Picks the measurement matching the user's chosen classification statistic.
function pick_stat(key, mean_val, min_val, max_val, median_val, stddev_val, intden_val, rawintden, area_px) {
    if (key == "mean") return mean_val;
    if (key == "median") return median_val;
    if (key == "min") return min_val;
    if (key == "max") return max_val;
    if (key == "stddev") return stddev_val;
    if (key == "integrated") return intden_val;
    if (key == "rawintegrated") return rawintden;
    if (key == "area") return area_px;
    return NaN;
}

function indexOfString(arr, val) {
    for (ii = 0; ii < arr.length; ii++) {
        if (arr[ii] == val) return ii;
    }
    return -1;
}


// Minimal CSV line parser that respects double-quoted fields (matches
// the simple quoting style written by the v1.5 segmentation pipeline).
// Returns an array of field strings for one line.
function parseCSVLine(line) {
    fields = newArray();
    nf = 0;
    cur = "";
    inQuotes = false;
    for (ci = 0; ci < lengthOf(line); ci++) {
        ch = substring(line, ci, ci + 1);
        if (ch == "\"") {
            inQuotes = !inQuotes;
        } else if (ch == "," && !inQuotes) {
            fields[nf] = cur;
            nf++;
            cur = "";
        } else {
            cur = cur + ch;
        }
    }
    fields[nf] = cur;
    nf++;
    return fields;
}
