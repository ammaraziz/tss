# script will ingest clades.tsv mutation define clades file (clades.tsv) defined by nextstrain, change the names to match internal gene naming
input_path = "../config/clades"
path = list.files(path = input_path, pattern = "*ha.tsv", full.names = T)

ha1_size = list('h1n1pdm' = 327,
                'h3n2' = 329,
                'vic' = 346)

convert = function(df, subtype){
  # change to HA numbering
  ha2_ind = which(df$gene %in% "HA2")
  df[ha2_ind, 'site'] =  df[ha2_ind, 'site'] + ha1_size[[subtype]]
  
  #rename
  df$gene = gsub(pattern = "HA\\d", replacement = "HA", df$gene)
  
  return(df)
}

get_subtype = function(file_path){
  gsub(pattern = "clades_(\\w+)_ha.tsv", replace = "\\1", basename(file_path))
}

for (f in path) {
  df = read.delim(f, sep = "\t")
  subtype = get_subtype(f)
  new = convert(df = df, subtype = subtype)
  write.table(x = new, 
              file = paste0(input_path, "/clades_", subtype, "_ha.mod.tsv"), 
              sep = "\t", 
              row.names = F)
}



