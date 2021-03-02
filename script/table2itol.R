#!/conda2/bin/Rscript --vanilla


################################################################################
#
# table2itol.R -- Rscript script for generating input files for iTOL.
#
# (C) since 2016 by Markus Goeker (markus [DOT] goeker [AT] dsmz [DOT] de)
#
# This script is distributed under the terms of the GNU General Public License.
# See http://www.gnu.org/licenses/gpl.html for further information.
#
# table2itol.R is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You must install R (https://cran.r-project.org/) and for non-interactive use
# also optparse (https://cran.r-project.org/package=optparse) to run this
# script. Several other R packages are needed for special purposes. See the
# README for details.
#
# This script was written for the command line but can also be used in
# interactive mode. See the README for details.
#
################################################################################


# Option processing
#
if (!interactive() || length(find.package("optparse", NULL, TRUE))) {

  optionparser <- optparse::OptionParser(option_list = list(

    optparse::make_option(opt_str = c("-a", "--abort"), action = "store_true",
      help = paste0("Abort if a requested column cannot be found instead of ",
        "just skipping the data set [default: %default]"),
      default = FALSE),

    optparse::make_option(opt_str = c("-b", "--background"), type = "character",
      help = paste0("Column to define the background colours of the tip ",
        "labels; empty means no background colours [default: %default]"),
      metavar = "NAME", default = ""),

    optparse::make_option(opt_str = c("-c", "--conversion"), type = "character",
      help = paste0("Convert integer columns to factors ('factor') or numbers ",
        "with decimal points ('double') or just not 0/1 to logical vectors ",
        "('keep') [default: %default]"),
      metavar = "NAME", default = "none"),

    optparse::make_option(opt_str = c("-C", "--colour-file"),
      help = paste0("File in YAML format defining alternative colour vectors ",
        "for domain output [default: %default]"),
      metavar = "FILE", default = "", type = "character"),

    optparse::make_option(opt_str = c("-d", "--double-to-bars"),
      help = paste0("Create bar charts, not gradients, from numbers ",
        "with decimal points ('double') [default: %default]"),
      action = "store_true", default = FALSE),

    optparse::make_option(opt_str = c("-D", "--directory"), type = "character",
      help = paste0("Place output files in this directory ('.' means working ",
        "directory, empty means input file directory) [default: %default]"),
      metavar = "DIR", default = "."),

    optparse::make_option(opt_str = c("-e", "--emblems"), type = "character",
      help = paste0("Column to define symbol assignments; ignored if empty ",
        "[default: %default]"),
      metavar = "NAME", default = ""),

    optparse::make_option(opt_str = c("-f", "--favour"), type = "numeric",
      help = paste0("Numeric factor for favouring colours over symbols ",
        "(higher => more colours relative to symbols) [default: %default]"),
      metavar = "NUMBER", default = 1),

    optparse::make_option(opt_str = c("-G", "--gradient-file"),
      help = paste0("File in YAML format defining alternative colours ",
        "for gradient and binary output [default: %default]"),
      type = "character", metavar = "FILE", default = ""),

    optparse::make_option(opt_str = c("-h", "--help"), action = "store_true",
      help = "Show this help message, then exit [default: %default]",
      default = FALSE),

    optparse::make_option(opt_str = c("-i", "--identifier"), type = "character",
      help = paste0("Mandatory identifier column; after modification ",
        "as defined by --template this column must yield the tip labels of ",
        "the tree [default: %default]"),
      metavar = "NAME", default = "ID"),

    optparse::make_option(opt_str = c("-j", "--identifier2"),
      help = paste0("Optional 2nd identifier column, causing output of branch",
        "symbols; together with -i this identifies a node [default: %default]"),
      metavar = "NAME", default = "", type = "character"),

    optparse::make_option(opt_str = c("-l", "--label"), type = "character",
      help = paste0("Column to define the tip labels displayed in the picture ",
        "in place of the tip labels found in the tree [default: %default]"),
      metavar = "NAME", default = "Label"),

    optparse::make_option(opt_str = c("-m", "--max-size"), type = "integer",
      help = paste0("Exceeding this threshold causes fewer colours and more ",
        "symbols to be selected (see also --favour); also determines size of ",
        "branch symbols [default: %default]"),
      metavar = "INTEGER", default = 20L),

    optparse::make_option(opt_str = c("-n", "--na-strings"), type = "character",
      help = paste0("Sentinels for missing input values; several can be ",
        "provided, separated by the value of --separator [default: %default]"),
      metavar = "TEXT", default = "\t(null)\tNA"),

    optparse::make_option(opt_str = c("-o", "--opacity"), type = "numeric",
      help = paste0("Numeric factor for the transparency of the colours ",
        "(0 => transparent, 1 => fully opaque) [default: %default]"),
      metavar = "NUMBER", default = 1),

    optparse::make_option(opt_str = c("-p", "--precision"), type = "integer",
      help = paste0("Number of decimal points used in the gradient legends ",
        "[default: %default]"),
      metavar = "INTEGER", default = 1L),

    optparse::make_option(c("-r", "--restrict"), type = "character",
      help = paste0("How to select from numeric values that yield branch ",
        "symbols [default: %default]"),
      metavar = "TEXT/NUMBER", default = ""),

    optparse::make_option(opt_str = c("-s", "--separator"), type = "character",
      help = paste0("Input column separator for CSV-like files ",
        "[default: %default]"),
      metavar = "CHARACTER", default = "\t"),

    optparse::make_option(opt_str = c("-t", "--template"), type = "character",
      help = paste0("Template for sprintf function to convert ID column when ",
        "deviating from tip labels [default: %default]"),
      metavar = "PATTERN", default = "%s"),

    optparse::make_option(opt_str = c("-w", "--width"), type = "numeric",
      help = paste0("Border with used for domains, colour strips etc. ",
        "[default: %default]"),
      metavar = "NUMBER", default = 0.5)

  ), add_help_option = FALSE, prog = "table2itol.R",
  usage = "%prog [options] file1 file2 ...", description = "
  %prog: converting spreadsheet files to iTOL input, version 2.5.1",
  epilogue = "
FREQUENTLY NEEDED OPTIONS:

-i\tUnless name of tip identifier column happens to match default.
-l\tUnless name of final tip label column happens to match default.
-s\tUnless separator character happens to match default.

USE OF DATA TYPES:

character, integer, logical -> factor -> iTOL domains
integer -> double -> iTOL gradient | iTOL branch symbols
integer[, double] -> iTOL simplebar
logical -> iTOL binary | iTOL collapsing instructions

EXAMPLES:

# Set identifier, label and label background column; also prepend
# 'T' to ID column 'Genome_ID' (which must contain integers):
'%prog -i Genome_ID -l Strain -b Phylum -t T%i annotation.tsv'

For more examples see the test folder and the FAQ.
"
  )

  invisible(list2env(optparse::parse_args(optionparser,
    commandArgs(TRUE), TRUE, TRUE), environment()))

  if (length(args) && !options$help) {
    options$infiles <- args
    options$help <- NULL
    names(options) <- chartr("-", ".", names(options))
  } else {
    optparse::print_help(optionparser)
    if (interactive()) {
      rm(optionparser, options, args)
      message("
********************************************************************************

Apparently this script is running in interactive mode. You could now generate
iTOL files by setting some 'infiles' variable to a vector of file names and then
calling:

create_itol_files(infiles)

********************************************************************************
      ")
    } else {
      quit("no", 1L)
    }
  }

}


################################################################################


# Main function, does everything given a vector of file names 'infiles'. The
# other arguments are optional.
#
create_itol_files <- function(infiles, identifier = "ID", label = "Label",
    background = "", identifier2 = "", directory = ".", colour.file = "",
    gradient.file = "", separator = "\t", na.strings = paste0(c("", "(null)",
      "NA"), collapse = separator), abort = FALSE, conversion = "none",
    double.to.bars = FALSE, emblems = "", template = "%s", max.size = 20L,
    favour = 1, width = 0.5, precision = 1L, restrict = "", opacity = 1) {


  OLDOPT <- options(warn = 1L)
  on.exit(options(OLDOPT))


  # EL  ellipse
  # RE  rectangle
  # TL  left pointing triangle
  # TR  right pointing triangle
  # DI  rhombus (diamond)
  # HH  horizontal hexagon
  # HV  vertical hexagon
  # PL  left pointing pentagram
  # PR  right pointing pentagram
  # PU  up pointing pentagram
  # PD  down pointing pentagram
  # OC  octagon
  # GP  rectangle (gap; black filled rectangle with 1/3 normal height)
  #
  SYMBOLS <- c("EL", "RE", "TL", "TR", "DI", "HH", "HV",
    "PL", "PR", "PU", "PD", "OC", "GP")


  # Restricted set: 1 = square, 2 = circle, 3 = star, 4 = right triangle,
  # 5 = left triangle.
  #
  BRANCH_SYMBOLS <- seq_len(5L)


  # Like branch symbols. We omit #6 (check mark) because it does not display
  # nicely.
  #
  BINARY_SYMBOLS <- seq_len(5L)


  BLACK <- "#000000"


  LIGHTGREY <- "#E5E5E5" # R's grey90


  WHITE <- "#FFFFFF"


  OUTPUT_SEPARATOR <- "\t"


  # Used as end points of colour gradients, with white as other end point,
  # and for colouring binary data (logical vectors)
  #
  SPECIAL_COLORS <- c("#1f78b4", "#e31a1c", "#33a02c", "#b15928",
    "#6a3d9a", "#ff7f00")


  # Colour vectors collected by Jan P. Meier-Kolthoff.
  #
  COLOURS <- list(
    # Dark2; colorblind-safe
    JMK01 = "#1b9e77",
    # Dark2; colorblind-safe
    JMK02 = c("#1b9e77", "#d95f02"),
    # Dark2; colorblind-safe
    JMK03 = c("#1b9e77", "#d95f02", "#7570b3"),
    # 4-class Paired; colorblind-safe
    JMK04 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
    # 5-class Accent; print-friendly
    JMK05 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99"),
    # 6-class Paired; print-friendly
    JMK06 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c"),
    # 7-class Paired; print-friendly
    JMK07 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c", "#fdbf6f"),
    # Dark2; print-friendly
    JMK08 = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"),
    # 9-class Set1; print-friendly
    JMK09 = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
      "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"),
    # 10-class Paired
    JMK10 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"),
    # 11-class Paired
    JMK11 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
      "#ffff99"),
    # 12-class Paired
    JMK12 = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
      "#ffff99", "#b15928"),
    ## from here on: iwanthue (all colors, hard)
    JMK13 = c("#8393c7", "#8ad256", "#6a49c5", "#d2b351",
      "#cb55c3", "#4d4040", "#c4527c", "#57743d", "#d85439", "#7accb1",
      "#925136", "#ceb2ab", "#512f67"),
    JMK14 = c("#a2d1cd", "#5d39a8", "#71d14c", "#cb56c7",
      "#7ed094", "#4d4040", "#7077b8", "#c28b4c", "#cd9dae", "#c64a34",
      "#55868c", "#cccb51", "#b2436e", "#567137"),
    JMK15 = c("#92d4ad", "#6842c1", "#6ecf58", "#cb4ec2",
      "#55733d", "#4d4040", "#c99447", "#9083cb", "#c9d14f", "#4d2c63",
      "#cea4a2", "#d54f38", "#71a6bd", "#ca507f", "#823f33"),
    JMK16 = c("#76a5bd", "#bfdf44", "#cf4bab", "#66c95b",
      "#7c42c5", "#4d4040", "#7279ca", "#c27837", "#4b2a62", "#c7b956",
      "#cc8cb5", "#536e3b", "#d74746", "#84d3ae", "#893b42", "#cdb19a"),
    JMK17 = c("#823f35", "#77d952", "#6d44c4", "#78d5a1",
      "#cf4a70", "#4d4040", "#ca53bd", "#69923c", "#6d7fc4", "#d1d04e",
      "#532b63", "#d64d31", "#4b623d", "#ca96b7", "#78b5c2", "#ccbf9b",
      "#c58741"),
    JMK18 = c("#697bc5", "#5e9742", "#6641c0", "#7bdc57",
      "#c954c9", "#4d4040", "#4d2b62", "#73d6ac", "#d6493d", "#75adbe",
      "#c54883", "#526339", "#caca9b", "#7b332e", "#cfcf49", "#c89dc8",
      "#c58738", "#c78980"),
    JMK19 = c("#9e693f", "#9147d5", "#c9d747", "#9482d3",
      "#61913d", "#4d4040", "#6dd85e", "#d049a4", "#76d0b6", "#d5493c",
      "#6897bb", "#d7993d", "#553291", "#c7cb8a", "#472f5b", "#cd7993",
      "#496340", "#ccb8bc", "#7f2c3a"),
    JMK20 = c("#7295c1", "#d44b38", "#6ad14f", "#6a3bc0",
      "#cedb44", "#4d4040", "#77d192", "#cb4fc3", "#b1b85f", "#7772cc",
      "#d9973b", "#4f2b62", "#79d1cf", "#cc497b", "#4a6c2e", "#c990b5",
      "#752e30", "#d1c5ac", "#a26f47", "#537e71"),
    JMK21 = c("#90b5d9", "#d6532d", "#c84ccc", "#74d147",
      "#512d79", "#4d4040", "#6740c8", "#cace49", "#6b79d1", "#6ccc84",
      "#c8478c", "#74c4b8", "#cc4458", "#4f6781", "#cb9142", "#552443",
      "#c6cb97", "#82442d", "#c489c5", "#546d37", "#cb9896"),
    JMK22 = c("#392c51", "#4d4040", "#642c79", "#792d3b",
      "#6a3ec6", "#875b30", "#4f7231", "#547f72", "#d24637", "#6d71ce",
      "#d2497e", "#cd4fc8", "#6a8fbc", "#d88742", "#c78dc6", "#cc9795",
      "#c7af40", "#68cd55", "#72d4a6", "#9ecfd6", "#c9cb8f", "#c3de48"),
    JMK23 = c("#8ad93f", "#c749c4", "#5e8f3d", "#6639be",
      "#73d979", "#4d4040", "#d4ca4a", "#6c6ccc", "#d78c3b", "#6485b9",
      "#d24635", "#70d4ae", "#cc4279", "#cbcb99", "#4c295f", "#ce867e",
      "#793130", "#84cbd7", "#896c35", "#c27bbb", "#364e27", "#cab2cb",
      "#5b837b"),
    JMK24 = c("#ccc79a", "#6a42c7", "#d0a540", "#cc49c9",
      "#6dd755", "#4d4040", "#de5a26", "#7cc7d0", "#cc3f47", "#78d8a5",
      "#5e2d78", "#c9da51", "#6679d0", "#bf7348", "#c6b7d8", "#5f903c",
      "#c47ec5", "#6a5b29", "#ce4684", "#497359", "#772d38", "#c3858c",
      "#352444", "#5b7a9e"),
    JMK25 = c("#6ba43c", "#c74ace", "#cbe14b", "#6847cd",
      "#6ede53", "#4d4040", "#cbb248", "#592e82", "#d6842f", "#5e78c1",
      "#76dd99", "#c6438e", "#4b8047", "#cf4c67", "#7acdc4", "#d2472f",
      "#7ba5c4", "#79322f", "#c388cf", "#78662f", "#45294d", "#c8cd9d",
      "#3e5d4a", "#d08c6c", "#c698a9"),
    JMK26 = c("#73d991", "#b44adb", "#71d94d", "#cf4cb4",
      "#ccde4d", "#4d4040", "#ceae44", "#5a41c2", "#cdd09c", "#652e7a",
      "#83d7ce", "#dc4338", "#536e83", "#d34a79", "#5d9073", "#c68dc7",
      "#619339", "#85b1d7", "#da8340", "#6978cb", "#9d4533", "#34284e",
      "#d09e9e", "#732d41", "#364e25", "#866a38"),
    JMK27 = c("#363258", "#6ed853", "#5b3fc7", "#c9de43",
      "#b54ad9", "#4d4040", "#5c2c7e", "#b7d17b", "#cf4a83", "#6ed9a4",
      "#cd4450", "#8fd3d5", "#d74527", "#769ac1", "#d27d3f", "#6d75cf",
      "#d4af42", "#4f8c3b", "#d14eba", "#568778", "#c692c8", "#344625",
      "#d4c7a6", "#722e4c", "#c88988", "#7a3a25", "#86783a"),
    JMK28 = c("#7f3a27", "#71da53", "#c14bd4", "#55933d",
      "#626ad0", "#4d4040", "#623ac4", "#cbd943", "#542c79", "#c1d483",
      "#bc7fd0", "#6ad7a3", "#d84330", "#71bec7", "#ce7537", "#6f99d8",
      "#d5aa43", "#546586", "#7c7233", "#ce429f", "#3e6344", "#ce7d9f",
      "#2d1d38", "#c6b3ce", "#793151", "#bfcbae", "#d24566", "#c8927d"),
    JMK29 = c("#cdc2c2", "#663dc8", "#76dd51", "#c64ece",
      "#cfda49", "#4d4040", "#549e3f", "#7577da", "#d3522e", "#7cd6ce",
      "#d4425b", "#77de9a", "#542a7e", "#d1d395", "#321e3d", "#d74a98",
      "#95963d", "#586095", "#db9a3e", "#77abd9", "#8b3c67", "#639575",
      "#d08982", "#456129", "#ca92cc", "#896134", "#597984", "#742c28",
      "#283a28"),
    JMK30 = c("#31223c", "#bbe141", "#c94edb", "#65d559",
      "#8b3899", "#4d4040", "#613ec8", "#df9b36", "#6e75d5", "#c16c39",
      "#402a74", "#cfc248", "#da47a4", "#63d6ad", "#d94330", "#6abccd",
      "#c58181", "#617fae", "#7f2f2c", "#b5cfb8", "#833b65", "#b5d888",
      "#cc88cb", "#4e8a3b", "#d6466a", "#476d58", "#d2b284", "#544320",
      "#c9b6d0", "#867c36"),
    JMK31 = c("#913d83", "#ced242", "#6643d0", "#79d949",
      "#c249d4", "#4d4040", "#db45a4", "#68dc88", "#3a1f4f", "#c3d483",
      "#532e8e", "#da983e", "#6d79d5", "#9b4b29", "#d085d5", "#8b7d3b",
      "#c9a0c0", "#54913d", "#dc4b32", "#72d4b1", "#8f3e58", "#90d0d8",
      "#592720", "#d2c7a9", "#21262c", "#d64769", "#3b4f25", "#6ea2cf",
      "#cd887a", "#5c6089", "#568477"),
    JMK32 = c("#8f8b38", "#663cc8", "#6bd546", "#c74cce",
      "#b1d773", "#4d4040", "#c6e03a", "#59287c", "#5edb86", "#d14592",
      "#7ad9b1", "#da4627", "#719cd8", "#dc973a", "#6e71d7", "#dbc348",
      "#ca84c8", "#4c8b3a", "#d5445a", "#84ccd6", "#7f3353", "#d3c99f",
      "#2e1c38", "#ca7442", "#5a558b", "#803325", "#537286", "#cc8585",
      "#314826", "#cab3cc", "#7e6136", "#618d75"),
    JMK33 = c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41",
      "#a152dd", "#4d4040", "#5139c2", "#ceaa3b", "#432d7c", "#c6d179",
      "#8f379a", "#70d68c", "#d9432f", "#6ad5be", "#d5416a", "#76c2d7",
      "#d87a71", "#6a75d5", "#836834", "#c988d1", "#598939", "#7a3260",
      "#bed3b3", "#8f372e", "#6082b3", "#d47c35", "#312749", "#d4ac8b",
      "#314825", "#cab9d7", "#4b211f", "#ad788b", "#568275"),
    JMK34 = c("#d8436c", "#653cc7", "#b4dc41", "#d143d0",
      "#5fd857", "#4d4040", "#a4db84", "#c64496", "#6adcad", "#de4830",
      "#6aa3d9", "#d98731", "#6271d1", "#dec841", "#b062cd", "#528e36",
      "#c28acd", "#675b2c", "#cbb7d3", "#a53332", "#528089", "#532878",
      "#d9d393", "#2a1e3c", "#8ed4d3", "#834629", "#5e5e8a", "#a08e3c",
      "#2b482a", "#d78763", "#619470", "#c87b8d", "#702944", "#c3a994"),
    JMK35 = c("#72d4cf", "#ccdf3e", "#5533c1", "#70d951",
      "#ac42d6", "#4d4040", "#6d66dc", "#b9c866", "#562a84", "#71da99",
      "#db43c7", "#518f39", "#d04497", "#314826", "#bc6cc9", "#5d8b74",
      "#d2416d", "#72abd3", "#dd461f", "#6078c6", "#d7ab3b", "#c49ad6",
      "#7d6b2f", "#cab8c4", "#3c1a20", "#c8ddb6", "#312652", "#cfb182",
      "#7c3463", "#c98271", "#576782", "#d24243", "#cb7a99", "#82372d",
      "#cf7734"),
    JMK36 = c("#6ade4b", "#6344d3", "#7bdc86", "#b746d4",
      "#65a234", "#4d4040", "#dbc941", "#552c93", "#bee148", "#dc3fb4",
      "#62d7b4", "#903a7e", "#4a8245", "#cf74d0", "#da993a", "#3e255f",
      "#c0d3b2", "#291d2d", "#cdce7e", "#752c41", "#7dcbd6", "#c43c44",
      "#669bcf", "#de4e28", "#5b5e83", "#c97449", "#bd92d0", "#847933",
      "#d7417a", "#558279", "#d07d92", "#364525", "#ceb9d0", "#763d23",
      "#6872d2", "#be9880"),
    JMK37 = c("#645b8e", "#80dc40", "#4f2ea4", "#69dc7b",
      "#d848cd", "#4d4040", "#8548da", "#c7d84e", "#96368e", "#afd995",
      "#d54227", "#61d9b9", "#db4187", "#4a9339", "#cd83d6", "#7a8431",
      "#6870d5", "#e3bc3b", "#6b9bd7", "#d87935", "#6fbfcf", "#cd3e50",
      "#c3d8c8", "#772e29", "#dbc38b", "#3f2267", "#bf9340", "#cab1d6",
      "#304726", "#b2918d", "#2a1f35", "#d5816f", "#5e8c6b", "#c77192",
      "#497080", "#7d592d", "#732d52"),
    JMK38 = c("#cf8ad0", "#74e042", "#b946da", "#5be080",
      "#5834c1", "#4d4040", "#d248bb", "#59a434", "#8064d4", "#b4dc4e",
      "#893876", "#96db99", "#d9478a", "#499052", "#627bcf", "#dfd238",
      "#47277a", "#908f39", "#79a2d8", "#d79234", "#4c7788", "#df502c",
      "#625984", "#d7d27b", "#2e1d3b", "#6bdac4", "#d34557", "#6a8b73",
      "#9e4427", "#cfb5cd", "#78562e", "#7cc6d5", "#26392b", "#cdcfb2",
      "#702735", "#bd7984", "#405924", "#d59571"),
    JMK39 = c("#8b308f", "#74dd41", "#6939ca", "#cce346",
      "#d545d2", "#4d4040", "#b271dd", "#e39b39", "#5050bc", "#cabc46",
      "#3a1f64", "#5cde7e", "#d9428e", "#57a56d", "#d63949", "#76dfc2",
      "#7e3052", "#b7e28f", "#d286c6", "#66a234", "#6d83d8", "#d65629",
      "#76c3d2", "#843326", "#6aa0d5", "#9c762c", "#5f5488", "#d48e70",
      "#4a6a81", "#d36778", "#466b2c", "#b28491", "#273825", "#c1b47a",
      "#301b31", "#d0d2bd", "#6c552d", "#c9b8d8", "#5f8675"),
    JMK40 = c("#3c2b5d", "#dee032", "#ab48d5", "#5bd749",
      "#db49c6", "#4d4040", "#5c42d0", "#a4e040", "#462687", "#d8b136",
      "#8d3989", "#60d076", "#d7468f", "#63d8b5", "#de4528", "#77c7d6",
      "#d13a55", "#5f8c7b", "#ce88d5", "#759b31", "#696ecd", "#de8739",
      "#6f9ad6", "#b75738", "#aadc90", "#946d89", "#d0dc6a", "#2c1a25",
      "#c6d8bc", "#782849", "#ceb977", "#283f27", "#d9798c", "#447c3d",
      "#ceb8d4", "#635b2d", "#c79783", "#733426", "#476682", "#98762e")
  )


  GENERATED_FILES <- new.env(TRUE, emptyenv())


  # Note that type.convert() does not recognize some spellings of false/true.
  # The value for TRUE must be the first one of each row.
  #
  BINARY_VALUE_SPELLINGS <- matrix(c(
     "y", "n",
     "t", "f",
     "true", "false",
     "yes", "no",
     "on", "off"
  ), 5L, 2L, TRUE)

  ## Helper functions


  # E.g. anyNA() is only available from 3.1.0 on.
  #
  check_R_version <- function(wanted = numeric_version("3.2.0")) {
    if (getRversion() < wanted)
      stop(sprintf("need a newer version of R, %s or higher", wanted))
    invisible(TRUE)
  }


  # Checking makes sense because colour vectors can be user-defined.
  #
  check_colour_vectors <- function(x, lengthcheck) {
    if (lengthcheck && !identical(seq_along(x), lengths(x, FALSE)))
      stop("incorrectly arranged colour vectors")
    bad <- !vapply(x, is.character, NA)
    if (any(bad))
      stop("wrong data type of colour vector(s) no. ",
        paste0(seq_along(x)[bad], collapse = ", "))
    bad <- vapply(x, anyDuplicated.default, 0L)
    if (any(bad))
      stop("duplicated value in colour vector(s) no. ",
        paste0(seq_along(x)[bad > 0L], collapse = ", "))
    invisible(TRUE)
  }


  standardize_colour <- function(x, opacity) {
    x <- grDevices::col2rgb(x, TRUE)
    x["alpha", ] <- as.integer(x["alpha", ] * opacity)
    tolower(grDevices::rgb(x["red", ], x["green", ], x["blue", ],
      if (all(x["alpha", ] == 255L)) NULL else x["alpha", ], NULL, 255L))
  }


  # For input of user-defined colour vectors.
  #
  read_colour_vectors <- function(file, upto) {
    if (!nzchar(file))
      return(NULL)
    x <- yaml::yaml.load_file(file)
    if (!is.list(x))
      x <- list(x)
    n <- lengths(x)
    x[n > 0L & n <= upto]
  }


  # Input method dispatch is based on file extension. Depends on extra library
  # for Excel and Libreoffice/Openoffice files, respectively.  Must ensure
  # character vectors are converted to factors.
  #
  read_file <- function(file, sep, na) {

    read_xl <- function(sheet, path, na) {
      # for some reason read_excel() yields a 'tibble' instead of a data frame,
      # which does not display the same behaviour of `[`; hence we convert
      tryCatch(expr = as.data.frame(readxl::read_excel(path = path, na = na,
        sheet = sheet, col_names = TRUE, col_types = NULL, skip = 0L,
        trim_ws = FALSE)), error = function(e) {
          warning(e) # a typical error is to encounter an empty sheet
          data.frame() # now we can treat this later on ourselves
        })
    }

    rescue_integers <- function(x) { # necessary for 'tibble' input
      is_whole_number <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
        all(is.na(x) | abs(x - round(x)) < tolerance) # see ?is.integer
      }
      for (i in which(vapply(x, is.double, NA)))
        if (is_whole_number(x[, i]))
          storage.mode(x[, i]) <- "integer"
      x
    }

    rescue_factors <- function(x) { # not necessary for CSV input
      for (i in which(vapply(x, is.character, NA)))
        x[, i] <- factor(x[, i])
      x
    }

    switch(
      EXPR = tolower(tools::file_ext(file)),
      ods = lapply(lapply(X = readODS::ods_sheets(file), path = file,
        FUN = readODS::read_ods, na = na[[1L]], col_names = TRUE,
        col_types = NULL, formula_as_formula = FALSE, skip = 0L, range = NULL),
        rescue_factors),
      xls =,
      xlsx = lapply(lapply(lapply(readxl::excel_sheets(file),
        read_xl, file, na), rescue_integers), rescue_factors),
      list(read.table(file = file, header = TRUE, sep = sep, quote = "\"",
        na.strings = na, fill = FALSE, stringsAsFactors = TRUE,
        dec = ".", check.names = FALSE, comment.char = ""))
    )

  }


  # Checking makes sense if 'chr' is used to join strings together.
  #
  assert_no_forbidden_character <- function(chr, ...) {
    x <- list(...)
    for (str in lapply(x[vapply(x, is.factor, NA)], levels.default)) {
      bad <- grepl(chr, str, FALSE, FALSE, TRUE)
      if (any(bad))
        stop(sprintf("string '%s' contains forbidden character '%s'",
          str[bad][[1L]], chr))
    }
    invisible(TRUE)
  }


  # We add white at the end, assuming this represents NA, when NA values occur.
  #
  select_colours <- function(size, hasna) {
    if (hasna) {
      message(sprintf("Fetching %i colour(s) ...", size - 1L))
      if (size < 2L)
        WHITE
      else
        c(COLOURS[[size - 1L]], WHITE)
    } else {
      message(sprintf("Fetching %i colour(s) ...", size))
      if (size < 1L)
        character()
      else
        COLOURS[[size]]
    }
  }


  # Used for generating legend titles.
  #
  pretty_str <- function(x) {
    chartr("_.", "  ", x)
  }


  # We assume NA values have already been removed.
  #
  legend_range <- function(x, precision) {
    if (length(precision))
      sprintf(sprintf("%%s (%%.%if)", precision), c("Min.", "Max."), range(x))
    else
      sprintf("%s (%i)", c("Min.", "Max."), range(x))
  }


  # Used to not display branch symbols associated with certain values.
  #
  mask_if_requested <- function(x, cutoff, restriction) {
    outliers <- function(x, n) {
      me <- median(x, na.rm = TRUE)
      ma <- mad(x, na.rm = TRUE)
      x > me + ma * n | x < me - ma * n
    }
    if (is.na(cutoff))
      return(logical(length(x)))
    switch(
      EXPR = restriction,
      atleast = x < cutoff,
      beyond = !outliers(x, cutoff),
      larger = x <= cutoff,
      smaller = x >= cutoff,
      upto = x > cutoff,
      within = outliers(x, cutoff),
      stop(sprintf("unkown 'restriction' value '%s'", restriction))
    )
  }


  # Used to modify several vectors at once. Needs at least one argument.
  #
  coordinated_na_removal <- function(...) {
    ok <- !is.na(..1)
    if (all(ok))
      return(FALSE)
    args <- list(...)
    names(args) <- all.names(match.call(), FALSE, -1L, FALSE)
    parentframe <- parent.frame()
    for (name in unique.default(names(args)))
      assign(name, args[[name]][ok], parentframe)
    TRUE
  }


  # Helper function to convert output vectors or matrices.
  #
  convert_to <- function(x, unavailable)  {
    storage.mode(x) <- typeof(unavailable)
    x[is.na(x)] <- unavailable
    x
  }


  # Used for generating the output filename.
  #
  itol_filename <- function(colname, kind, directory) {
    result <- file.path(directory, sprintf("iTOL_%s-%s.txt", kind,
      gsub("\\W", "_", colname, FALSE, TRUE)))
    if (exists(result, GENERATED_FILES))
      stop(sprintf("name clash: file '%s' has already been generated", result))
    GENERATED_FILES[[result]] <- TRUE
    message(sprintf("Generating %s file for column '%s' ...", kind, colname))
    result
  }


  # Here '...' contains the data part.
  #
  print_itol <- function(outdir, title, annotation, ...) {

    join <- function(x) {
      if (!length(x))
        return(NULL)
      if (is.null(names(x)))
        stop("non-empty annotation lists must have names")
      if (!all(vapply(x, is.atomic, NA)))
        stop("non-empty annotation lists must contain only atomic values")
      x <- x[sort.list(names(x))]
      sizes <- lengths(x, FALSE)
      for (i in which(sizes > 1L))
        x[[i]] <- paste0(x[[i]], collapse = OUTPUT_SEPARATOR)
      for (i in which(!sizes))
        x[[i]] <- ""
      paste(names(x), unlist(x, FALSE, FALSE), sep = OUTPUT_SEPARATOR)
    }

    if (is.character(annotation)) {
      colname <- annotation
      annotation <- NULL
    } else {
      colname <- get("DATASET_LABEL", annotation)
    }

    kind <- switch(
      EXPR = title,
      branchsymbols = "DATASET_SYMBOL",
      collapse = "COLLAPSE",
      labels = "LABELS",
      treecolors = "TREE_COLORS",
      binary =,
      colorstrip =,
      domains =,
      gradient =,
      heatmap =,
      simplebar =,
      text = sprintf("DATASET_%s", toupper(title)),
      stop(sprintf("unknown title '%s'", title))
    )

    separator <- sprintf("SEPARATOR %s", switch(
      EXPR = OUTPUT_SEPARATOR,
      `\t` = "TAB",
      stop(sprintf("output separator '%s' not yet supported", OUTPUT_SEPARATOR))
    ))

    file <- itol_filename(colname, title, outdir)

    cat(kind, separator, join(annotation), "DATA", file = file,
      labels = NULL, sep = "\n", fill = FALSE, append = FALSE)
    cat(paste(..., sep = OUTPUT_SEPARATOR, collapse = NULL), file = file,
      labels = NULL, sep = "\n", fill = FALSE, append = TRUE)

  }


  ## Functions for special columns


  # For labelling the leaves.
  #
  emit_itol_labeltexts <- function(x, ids, name, outdir, ...) {
    coordinated_na_removal(x, ids)
    print_itol(outdir, "labels", name, ids, x)
  }


  # For colouring the leaves. 'x' is a factor or coerced to a factor, hence NAs
  # do not get removed.
  #
  emit_itol_labelcolors <- function(x, ids, name, outdir, ...) {
    if (!is.factor(x))
      x <- factor(x)
    x <- addNA(x, TRUE)
    size <- length(levels.default(x))
    if (size > length(COLOURS)) {
      warning(sprintf("skipping column '%s', which yields > %i levels",
        name, length(COLOURS)))
      return()
    }
    annotation <- list(
      COLOR = "#a6cee3",
      DATASET_LABEL = name,
      LEGEND_COLORS = select_colours(size, anyNA(levels.default(x))),
      LEGEND_LABELS = levels.default(x),
      LEGEND_SHAPES = rep.int(1L, size),
      LEGEND_TITLE = pretty_str(name)
    )
    print_itol(outdir, "treecolors", annotation,
      ids, "range", annotation$LEGEND_COLORS[x], x)
  }


  ## Functions for columns according to data type (class)


  # Output varies depending on the number of colours and symbols chosen and/or
  # available. 'x' is a factor, hence NAs do not get removed.
  #
  emit_itol_factor <- function(x, ids, name, outdir, symbols, maxclrs,
      favour, borwid, ...) {

    product <- function(x, y) {
      cbind(rep(x = x, each = length(y)), rep.int(y, length(x)))
    }

    x <- addNA(x, TRUE)
    annot1 <- list(DATASET_LABEL = name, MARGIN = 5, COLOR = "#bebada")
    size <- length(levels.default(x))

    if (size > maxclrs * length(SYMBOLS)) {

      # additional columns: position, color, style, size_factor, rotation
      print_itol(outdir, "text", annot1, ids, x,
        -1, BLACK, "normal", 0.75, 0)

    } else if (length(symbols) || size > maxclrs) {

      if (length(symbols)) {
        symbols <- vapply(split.default(symbols, x), `[[`, "", 1L)
        clrs <- select_colours(size, anyNA(levels.default(x)))
      } else {
        nsym <- ncls <- ceiling(sqrt(size))
        nsym <- round(nsym / favour, 0L)
        ncls <- round(ncls * favour, 0L)
        if (nsym > length(SYMBOLS) || ncls > maxclrs) {
          msg <- sprintf(
            "Column '%s': # symbols (%i) or # colours (%i) inacceptable",
            name, nsym, ncls)
          if (favour >= 1) {
            ncls <- maxclrs
            nsym <- ceiling(size / ncls)
          } else {
            nsym <- length(SYMBOLS)
            ncls <- ceiling(size / nsym)
          }
          message(msg, sprintf(", trying %i/%i instead.", nsym, ncls))
        }
        clrs <- select_colours(ncls, FALSE) # NA treated below
        symbols <- product(SYMBOLS, clrs)
        clrs <- symbols[, 2L]
        if (anyNA(levels.default(x))) # ensure last selected position is white
          clrs[[size]] <- WHITE
        symbols <- symbols[, 1L]
      }

      annotation <- c(annot1, list(
        BACKBONE_HEIGHT = 0, # controls the height of the midline
        BACKBONE_COLOR = WHITE, # controls the color of the midline
        # we are hiding it by drawing it white
        BORDER_WIDTH = borwid,
        HEIGHT_FACTOR = 1,
        LEGEND_COLORS = clrs[seq_len(size)],
        LEGEND_LABELS = levels.default(x),
        LEGEND_SHAPES = symbols[seq_len(size)],
        LEGEND_TITLE = pretty_str(name),
        SHOW_DOMAIN_LABELS = 0,
        WIDTH = 25
      ))
      assert_no_forbidden_character("|", x)
      joint <- paste(symbols[x], 0L, 10L, clrs[x], as.character(x), sep = "|")
      print_itol(outdir, "domains", annotation, ids, 10L, joint)

    } else {

      annotation <- c(annot1, list(
        BORDER_WIDTH = borwid,
        LEGEND_COLORS = select_colours(size, anyNA(levels.default(x))),
        LEGEND_LABELS = levels.default(x),
        LEGEND_SHAPES = rep.int(1L, size),
        LEGEND_TITLE = pretty_str(name),
        STRIP_WIDTH = 25
      ))
      print_itol(outdir, "colorstrip", annotation,
        ids, annotation$LEGEND_COLORS[x], x)

    }
  }


  # Integer vectors yield a bar chart.
  #
  emit_itol_integer <- function(x, ids, name, outdir, precision, ...) {
    if (!is.object(x)) # then we can assume it is really an integer vector
      precision <- NULL
    coordinated_na_removal(x, ids)
    annotation <- list(
      COLOR = BLACK,
      DATASET_LABEL = name,
      LEGEND_COLORS = BLACK,
      LEGEND_LABELS = paste0(legend_range(x, precision), collapse = ", "),
      LEGEND_SHAPES = 1L,
      LEGEND_TITLE = pretty_str(name),
      MARGIN = 5,
      WIDTH = 200
    )
    print_itol(outdir, "simplebar", annotation, ids, x)
  }


  # For treating vectors of mode 'double' like integers.
  #
  emit_itol_pseudointeger <- function(...) {
    emit_itol_integer(...)
  }


  # Should not normally occur in input.
  #
  emit_itol_list <- function(x, ids, name, outdir, ...) {
    message(sprintf("Skipping column '%s' of mode 'list' ...", name))
  }


  # For logical vectors. NAs do not get removed.
  #
  emit_itol_logical <- function(x, ids, name, outdir, bincolor, binsymbol,
      borwid, ...) {
    annotation <- list(
      BORDER_WIDTH = borwid,
      COLOR = "#4daf4a",
      DATASET_LABEL = name,
      FIELD_COLORS = bincolor,
      FIELD_LABELS = pretty_str(name),
      FIELD_SHAPES = binsymbol,
      LEGEND_COLORS = bincolor,
      LEGEND_LABELS = pretty_str(name),
      LEGEND_SHAPES = 1L,
      LEGEND_TITLE = pretty_str(name),
      MARGIN = 5,
      WIDTH = 20
    )
    print_itol(outdir, "binary", annotation, ids, convert_to(x, -1L))
  }


  # Vectors of mode 'double' (of class 'numeric' in R) yield a colour gradient.
  #
  emit_itol_numeric <- function(x, ids, name, outdir, endcolor,
      precision, borwid, ...) {
    coordinated_na_removal(x, ids)
    annotation <- list(
      BORDER_WIDTH = borwid,
      COLOR = "#fb9a99",
      COLOR_MAX = endcolor,
      COLOR_MIN = LIGHTGREY,
      DATASET_LABEL = name,
      LEGEND_COLORS = c(LIGHTGREY, endcolor),
      LEGEND_LABELS = legend_range(x, precision),
      LEGEND_SHAPES = c(1L, 1L),
      LEGEND_TITLE = pretty_str(name),
      MARGIN = 5,
      STRIP_WIDTH = 50
    )
    print_itol(outdir, "gradient", annotation, ids, x)
  }


  emit_itol_matrix <- function(x, ids, name, outdir, endcolor,
      precision, borwid, ...) {
    if (is.integer(x))
      precision <- 0L
    annotation <- list(
      BORDER_WIDTH = borwid,
      COLOR = "#fb9a99",
      COLOR_MAX = endcolor,
      COLOR_MIN = LIGHTGREY,
      COLOR_NAN = WHITE,
      DATASET_LABEL = name,
      FIELD_LABELS = colnames(x),
      LEGEND_COLORS = c(LIGHTGREY, endcolor),
      LEGEND_LABELS = legend_range(x, precision),
      LEGEND_SHAPES = c(1L, 1L),
      LEGEND_TITLE = pretty_str(name),
      MARGIN = 5,
      STRIP_WIDTH = 50
    )
    x <- apply(X = convert_to(x, "X"), MARGIN = 1L,
      FUN = paste0, collapse = OUTPUT_SEPARATOR)
    print_itol(outdir, "heatmap", annotation, ids, x)
  }


  # Not yet implemented.
  #
  emit_branch_annotation_factor <- function(x, ids, name, outdir, ...) {
    message(sprintf("Skipping column '%s' of mode 'factor' ...", name))
  }


  # Not yet implemented.
  #
  emit_branch_annotation_integer <- function(x, ids, name, outdir, ...) {
    message(sprintf("Skipping column '%s' of mode 'integer' ...", name))
  }


  # Should not occur in input.
  #
  emit_branch_annotation_list <- function(x, ids, name, outdir, ...) {
    message(sprintf("Skipping column '%s' of mode 'list' ...", name))
  }


  # Vectors of class 'logical' select subtrees to be collapsed.
  #
  emit_branch_annotation_logical <- function(x, ids, name, outdir, ...) {
    x[is.na(x)] <- FALSE
    print_itol(outdir, "collapse", name, ids[x])
  }


  # Vectors of mode 'double' (of class 'numeric' in R) yield a colour gradient
  # within the branch symbols.
  #
  emit_branch_annotation_numeric <- function(x, ids, name, outdir, symbol,
      endcolor, branchpos, maxsize, precision, cutoff, restriction, ...) {
    coordinated_na_removal(x, ids)
    annotation <- list(
      COLOR = endcolor,
      DATASET_LABEL = name,
      LEGEND_TITLE = pretty_str(name),
      LEGEND_SHAPES = c(symbol, symbol),
      LEGEND_COLORS = c(LIGHTGREY, endcolor),
      LEGEND_LABELS = legend_range(x, precision),
      MAXIMUM_SIZE = maxsize
    )
    xclrs <- plotrix::color.scale(x = x, extremes = annotation$LEGEND_COLORS)
    mask <- mask_if_requested(x, cutoff, restriction)
    if (any(mask)) {
      x[mask] <- NA_real_
      coordinated_na_removal(x, xclrs, ids)
    }
    print_itol(outdir, "branchsymbols", annotation,
      # columns: ID, symbol, size, colour, fill, position
      ids, symbol, maxsize, xclrs, 1L, branchpos)
  }


  # Useful when R does not get the type right because of special notations;
  # also for user-defined type modifications.
  #
  fix_column_types <- function(x, convint, convdbl) {

    # convert binary integer vectors to logical vectors
    if (!is.element(convint, c("keep", "double")))
      for (i in which(vapply(x, is.integer, NA)))
        if (all(is.element(x[, i], c(0L, 1L, NA_integer_))))
          storage.mode(x[, i]) <- "logical"

    # convert factors to logical vectors if values look like boolean values
    if (convint != "keep")
      for (i in which(vapply(x, is.factor, NA))) {
        values <- tolower(levels.default(x[, i]))
        truevalue <- NA_character_
        for (j in seq_len(nrow(BINARY_VALUE_SPELLINGS)))
          if (all(is.element(values, BINARY_VALUE_SPELLINGS[j, ]))) {
            truevalue <- BINARY_VALUE_SPELLINGS[j, 1L]
            break
          }
        if (!is.na(truevalue))
          x[, i] <- tolower(x[, i]) == truevalue
      }

    if (convdbl)
      for (i in which(vapply(x, is.double, NA)))
        class(x[, i]) <- "pseudointeger"

    # convert integers and logical vectors to other data types if requested
    switch(
      EXPR = convint,
      keep = NULL,
      none = for (i in which(vapply(x, is.logical, NA)))
        if (anyNA(x[, i]))
          x[, i] <- factor(x[, i]),
      factor = {
        for (i in which(vapply(x, is.integer, NA)))
          x[, i] <- factor(x[, i])
        for (i in which(vapply(x, is.logical, NA)))
          x[, i] <- factor(x[, i])
      },
      double = {
        for (i in which(vapply(x, is.integer, NA)))
          storage.mode(x[, i]) <- "double"
        for (i in which(vapply(x, is.logical, NA)))
          x[is.na(x[, i]), i] <- FALSE
      },
      stop(sprintf("invalid integer/logical vector conversion indicator '%s'",
        convint))
    )

    x
  }


  # Helper function for itol_files().
  #
  assort <- function(f, x) {
    idx <- split.default(seq_along(f), f)
    if (length(x)) {
      result <- vector(typeof(x), length(f))
      for (i in idx)
        result[i] <- rep_len(x, length(i))
    } else {
      result <- vector("double", length(f))
      for (i in idx)
        result[i] <- seq_along(i) / (length(i) + 1L)
    }
    result
  }


  # Helper function for itol_files().
  #
  get_col <- function(name, x, strict) {
    if (length(name) != 1L)
      stop("need a single column name for identifying special column")
    result <- match(name, names(x), 0L)
    if (!result)
      if (strict)
        stop(sprintf(
          "selected column '%s' does not exist -- must select one of %s",
          name, paste0(sprintf("'%s'", names(x)), collapse = ", ")))
      else
        warning(sprintf("cannot find column '%s', skipping data", name))
    result
  }


  # Called by itol_files() instead of jumping to the main branch.
  #
  branch_annotation_files <- function(x, icol, jcol, idpat, precision, outdir,
      maxsize, restrict) {

    if (nzchar(restrict)) {
      cutoff <- as.double(basename(restrict))
      restriction <- dirname(restrict)
      if (restriction == ".") # when only a number was provided
        restriction <- "upto"
    } else {
      cutoff <- NA_real_
      restriction <- "upto"
    }

    idpos <- get_col(icol, x, TRUE)
    jpos <- get_col(jcol, x, TRUE)
    assert_no_forbidden_character("|", x[, idpos], x[, jpos])
    icol <- ifelse(is.na(x[, jpos]), sprintf(idpat, x[, idpos]),
      paste(sprintf(idpat, x[, idpos]), sprintf(idpat, x[, jpos]), sep = "|"))
    x <- x[, -c(idpos, jpos), drop = FALSE]

    klass <- vapply(x, class, "")

    # normal columns, dispatch done according to data type (class)
    mapply(FUN = function(fun, ...) fun(...), x = x, name = names(x),
      fun = lapply(sprintf("emit_branch_annotation_%s", klass), match.fun),
      endcolor = assort(klass, SPECIAL_COLORS),
      branchpos = assort(klass, NULL), SIMPLIFY = FALSE,
      symbol = assort(klass, BRANCH_SYMBOLS), USE.NAMES = FALSE,
      MoreArgs = list(ids = icol, precision = precision, outdir = outdir,
        maxsize = maxsize, cutoff = cutoff, restriction = restriction))

    invisible(TRUE)
  }


  ## Main


  # The main function, taking care of all columns of data frame 'x'.
  #
  itol_files <- function(x, lcol, bcol, icol, jcol, scol, idpat, precision,
      maxsize, favour, strict, convint, convdbl, outdir, borwid, restrict) {

    # identifier column (mandatory in strict mode), step 1
    idpos <- get_col(icol, x, strict)
    if (!idpos)
      return(invisible(FALSE))
    names(idpos) <- names(x)[[idpos]]

    if (anyNA(x[, idpos]))
      x <- x[!is.na(x[, idpos]), , drop = FALSE]

    if (!all(dim(x))) {
      if (strict)
        stop("encountered empty data frame")
      else
        warning("skipping empty data frame")
      return(invisible(FALSE))
    }

    x <- fix_column_types(x, convint, convdbl)

    # must be done before the first use of 'outdir'
    if (!dir.exists(outdir))
      dir.create(outdir)

    # generate branch symbols, skip normal run
    if (length(jcol) && all(nzchar(jcol)))
      return(branch_annotation_files(x = x, icol = icol, jcol = jcol,
        idpat = idpat, precision = precision, outdir = outdir,
        maxsize = maxsize, restrict = restrict))

    # identifier column (mandatory in strict mode), step 2
    icol <- x[, idpos]
    if (is.factor(icol))
      icol <- as.character(icol)
    icol <- sprintf(idpat, icol)

    # label column (mandatory in strict mode)
    lpos <- get_col(lcol, x, strict)
    if (lpos)
      emit_itol_labeltexts(x = x[, lpos], ids = icol,
        name = names(x)[[lpos]], outdir = outdir)

    # background colour column (optional in strict mode)
    if (length(bcol) && all(nzchar(bcol))) {
      cpos <- get_col(bcol, x, strict)
      if (cpos)
        emit_itol_labelcolors(x = x[, cpos], ids = icol,
          name = names(x)[[cpos]], outdir = outdir)
    } else {
      cpos <- 0L
    }

    # symbol-defining column (optional in strict mode)
    symbols <- NULL
    if (length(scol) && all(nzchar(scol))) {
      spos <- get_col(scol, x, strict)
      if (spos) {
        symbols <- x[, spos]
        if (!is.factor(symbols) || anyNA(symbols) ||
            length(levels.default(symbols)) > length(SYMBOLS)) {
          warning("column '", scol, "' is either not a factor or ",
            "has too many levels to be used for deriving symbols")
          symbols <- NULL
        } else {
          symbols <- SYMBOLS[symbols]
        }
      }
    }

    # normal columns, dispatch done according to data type (class)
    x <- x[, -c(idpos, lpos, cpos), drop = FALSE]
    klass <- vapply(x, class, "")

    # join all remaining columns in a matrix when applicable
    if (all(duplicated.default(klass)[-1L]))
      switch(
        EXPR = klass[[1L]],
        integer =,
        numeric = { # 'pseudointeger' should not be accepted here
          klass <- "matrix"
          x <- list(as.matrix(x))
          names(x) <- names(idpos)
        },
        NULL
      )

    clrs <- assort(klass, SPECIAL_COLORS)
    mapply(FUN = function(fun, ...) fun(...), x = x, name = names(x),
      fun = lapply(sprintf("emit_itol_%s", klass), match.fun), endcolor = clrs,
      bincolor = clrs, binsymbol = assort(klass, seq_along(BINARY_SYMBOLS)),
      MoreArgs = list(ids = icol, precision = precision, outdir = outdir,
        symbols = symbols, maxclrs = maxsize, favour = favour,
        borwid = borwid), SIMPLIFY = FALSE, USE.NAMES = FALSE)

    invisible(TRUE)

  }

  check_R_version()

  # assignment of input colour vectors is solely by vector length
  for (clrs in read_colour_vectors(colour.file, length(COLOURS)))
    COLOURS[[length(clrs)]] <- clrs
  COLOURS[] <- lapply(COLOURS, standardize_colour, opacity)
  check_colour_vectors(COLOURS, TRUE)

  # any length allowed, last one wins
  for (clrs in read_colour_vectors(gradient.file, Inf))
    SPECIAL_COLORS <- clrs
  SPECIAL_COLORS <- standardize_colour(SPECIAL_COLORS, opacity)
  check_colour_vectors(list(SPECIAL_COLORS), FALSE)

  LIGHTGREY <- standardize_colour(LIGHTGREY, opacity)

  na.strings <- unlist(strsplit(na.strings, separator, TRUE), FALSE, FALSE)
  if (!length(na.strings))
    na.strings <- ""

  for (infile in infiles)
    # note that read_file() is supposed to return a list of data frames
    lapply(X = read_file(infile, separator, na.strings), FUN = itol_files,
      bcol = background, precision = precision, lcol = label, icol = identifier,
      scol = emblems, idpat = template, maxsize = max.size, favour = favour,
      outdir = if (nzchar(directory)) directory else dirname(infile),
      strict = abort, jcol = identifier2, borwid = width, restrict = restrict,
      convint = conversion, convdbl = double.to.bars)

  invisible(NULL)

}


if (!interactive()) {
  do.call(create_itol_files, options)
}


################################################################################

