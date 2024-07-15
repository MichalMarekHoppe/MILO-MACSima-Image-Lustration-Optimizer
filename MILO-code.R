## MILO - MACSima Image Lustration Optimizer version 1.2.0
## by Dr Michal Marek Hoppe, 2024
## version 1.2.0

# Execute the code to organize preprocessed Miltenyi MACSima hyperplex-IHC images 

# currently supports only one-slide + one-well output
# i.e. files from a single project directory wuth R1 and A1 sub-directories

# tune parameters
vars <-
  list(comp = c("dapi", "mono"), # exports DAPI composite and monchrome images
       norm = c("nrm", "raw"), # exports normalized and absolute intensities
       jpg = c(15, 90), # exports compressed and relatively uncompressed jpegs
       format = c("ROIs", "markers")) # creates file names for ROI scrolling and marker scrolling

if (T) {
  invisible(library(magick))
  rm(list = ls()); gc()
  
  # functions
  z <- 
    function(x, d, what = " ", side = "right") {
      if (side == "right") {
        output <-
          paste0(x,
                 paste(rep(what, length.out = d - nchar(x)), collapse = ""))
      } else if (side == "left") {
        output <-
          paste0(paste(rep(what, length.out = d - nchar(x)), collapse = ""),
                 x)
      }
      return(output)
    }
  project <-
    gsub("./", "", list.dirs())[2]
  
  # select folder
  folder <- 
    list.dirs(paste0(project, "/"))[grep("ROI", list.dirs(paste0(project, "/")), invert = TRUE)]
  folder <- 
    paste0(folder[nchar(folder) == max(nchar(folder))], "/")
  
  rois <-
    gsub(folder, "", list.dirs(folder))[gsub(folder, "", list.dirs(folder)) != ""]
  
  # extract markers
  roi <- rois[1]
  markers <- 
    list.files(paste0(folder, "/", roi))[grep(".tif", list.files(paste0(folder, "/", roi)))]
  markers <-
    gsub(".tif", "", gsub(".*_A-", "",  gsub("_C-.*", "", markers)))
  markers <-
    unique(markers[!markers %in% c("None", "dummy")])
  names(markers) <-
    c(rep("00", length.out = grep("DAPI", markers)),
      formatC(1:(length(markers) - grep("DAPI", markers)),
              flag = "0", width = 2))
  fluors <-
    c("APC", "FITC", "PE")
  
  # create milo directory
  dirs <- 
    apply(as.matrix(expand.grid(paste0("milo-", project), "/", 
                                vars$norm, "/",
                                vars$comp, "/",
                                paste0("jpg", vars$jpg), "/",
                                vars$format)),
          1,
          paste, collapse = "")
  
  for (d in dirs) {
    dir.create(d,
               recursive = TRUE, 
               showWarnings = FALSE); 
  }
  
  # gear for website
  dir.create("website/img", 
             recursive = TRUE,
             showWarnings = FALSE)
  write.table(sort(ifelse(markers %in% fluors,
                          paste0("af-", markers),
                          markers))[c(grep(paste(fluors, 
                                                 collapse = "|"),
                                           markers, invert = TRUE),
                                      grep(paste(fluors, 
                                                 collapse = "|"),
                                           markers, invert = FALSE))],
              file = paste0(getwd(), "/website/markers.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(sort(rois),
              file = paste0(getwd(), "/website/ids.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # deal with images
  variations <- unlist(lapply(vars, length))
  variations[names(variations) == "format"] <-
    variations[names(variations) == "format"] + 1 # account for website export
  total <-
    length(rois) * length(markers) * prod(variations)
}
# define image size manually
width <- 6363; height <- 5286; x <- 0; r <- rois[1]
for (r in rev(rois)) {
  files <-
    list.files(paste0(folder, r, "/"))[grep(".tif",
                                            list.files(paste0(folder, r, "/")))]
  
  file_DAPI <-
    files[grep("^DAPI$", gsub(".tif", "", gsub(".*_A-", "",  gsub("_C-.*", "", files))))]
  if (length(file_DAPI) > 1) {
    file_DAPI <- file_DAPI[1]
  }
  #load DAPI
  img_dapi <-
    image_composite(image_convert(image_normalize(image_read(paste0(folder, r, "/", file_DAPI))),
                                  colorspace = "sRGB"),
                    image_blank(width = width, height = height, color = "#4F67AB"),
                    operator = "multiply")
  cat(paste0("Loaded DAPI: ", file_DAPI), "\n")
  m <- "Ki67"
  for (m in rev(markers)) {
    gc()
    file <-
      files[grep(paste0("^",m, "$"), gsub(".tif", "", gsub(".*_A-", "",  gsub("_C-.*", "", files))))]
    if (length(file) > 1) {
      file <- file[1]
    }
    # load marker image
    img_raw <-
      image_read(paste0(folder, r, "/", file))
    
    img_nrm <-
      image_normalize(img_raw)
    # merge with DAPI
    img_raw_merge <-
      image_composite(image_convert(img_raw, colorspace = "sRGB"),
                      img_dapi,
                      operator = "lighten")
    
    img_nrm_merge <-
      image_composite(image_convert(img_nrm, colorspace = "sRGB"),
                      img_dapi,
                      operator = "lighten")
    cat(paste0("Prepared images from: ", file), "\n")
    c <- vars$comp[1]
    for (c in vars$comp) {
      n <- vars$norm[1]
      for (n in vars$norm) {
        if        (c == "dapi" & n == "nrm") {
          export <- img_nrm_merge
        } else if (c == "dapi" & n == "raw") {
          export <- img_raw_merge
        } else if (c == "mono" & n == "nrm") {
          export <- img_nrm
        } else if (c == "mono" & n == "raw") {
          export <- img_raw
        }
        j <- vars$jpg[1]
        for (j in vars$jpg) {
          # website first
          filename <-
            paste0("website/img/",
                   paste0(gsub("[^A-Za-z]", "", r), # name for ROI focus
                          z(gsub("[^0-9]", "", r), 2, side = "left", what = 0), "_",
                          ifelse(m %in% c("APC", "FITC", "PE"),
                                 paste0("af-", m),
                                 m)), "_",
                   toupper(n), "_", toupper(c), "_",
                   "q", formatC(j, width = 3, flag = 0),
                   ".jpeg")
          image_write(export,
                      filename,
                      format = "jpeg", quality = j)
          gc()
          cat(paste0("Saved: ", filename), "\n")
          x <- x + 1
          f <- vars$format[1]
          for (f in vars$format) {
            file.copy(from = filename,
                      to = paste0(paste0("milo-", project), "/", n, "/", c, "/",
                                  "jpg", j, "/",
                                  f, "/",
                                  ifelse(f == "markers",
                                         paste0(names(markers)[grep(paste0("^",m, "$"),
                                                                    markers)], "_", # name for marker focus
                                                ifelse(m %in% c("APC", "FITC", "PE"),
                                                       paste0("af-", m),
                                                       m), "_",
                                                gsub("[^A-Za-z]", "", r),
                                                z(gsub("[^0-9]", "", r), 2, side = "left", what = 0)),
                                         paste0(gsub("[^A-Za-z]", "", r), # name for ROI focus
                                                z(gsub("[^0-9]", "", r), 2, side = "left", what = 0), "_",
                                                names(markers)[grep(paste0("^",m, "$"),
                                                                    markers)], "_",
                                                ifelse(m %in% c("APC", "FITC", "PE"),
                                                       paste0("af-", m),
                                                       m))),
                                  ".jpeg"),
                      overwrite = TRUE)
            x <- x + 1
            cat(paste0("Saved on ", z(r, 5), " ", z(m, max(nchar(markers))), " ",
                       z(x, nchar(total), side = "left"), "/", total, " ",
                       formatC(round(x/total * 100, 2), format = "f", digits = 2), "%"), "\r")
          }
        }
      }
    }
  }
}
