args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("No arguments!", call.=FALSE)
}

db_LOLA = args[1]
regions_of_interest = args[2]
bg_regions = args[3]
output_file = args[4]

message("DB: ", db_LOLA, " Input: ", regions_of_interest, " BG: ", bg_regions, " Output: ", output_file)

library("LOLA")
regionDB = loadRegionDB(db_LOLA)
regions = readBed(regions_of_interest)
bg = readBed(bg_regions)

checkUniverseAppropriateness(regions, bg, cores = 1, fast = FALSE)
locResults = runLOLA(regions, bg, regionDB, cores=1)

write.csv(locResults, output_file)
message("DONE")
