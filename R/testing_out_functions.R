###############################################################################
## Fiddling with the blast results for the Abars ##############################
###############################################################################


require(dplyr, quietly = TRUE)


## Functions 

get_input <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  if(length(args) != 2){
      print("Not correct number of inputs, need 2: <right_csv> <contig_loc>")
      quit(save = "no", status = 1, runLast = FALSE)
    }
    
  left_csv <- args[1]
  contig_dir <- args[2]
  
  args <- c(left_csv, contig_dir)
  return(args)
}
  
contig_location_adder <- function(isolate, contig_dir, hit_table, iter){
  ## Function to take an input isolate and add in the contig 
  ## numbers for the different hit locations 
  #if(iter == 71) browser()
  contig_dir <- gsub("/$","",contig_dir)
  contig_file <- paste(contig_dir , "/" , isolate , "#contig_bounds.csv", sep = "")

  contig_table <- read.csv(contig_file, stringsAsFactors = FALSE)
  hit_table$contig <- NA
  for(row in 1:nrow(hit_table)){
    current_hit <- hit_table[row,]
    current_pos <- c(current_hit$sstart, current_hit$send)
    contig <- which((contig_table[,1] <= (min(current_pos) + 15)) & ((contig_table[,2] >= (max(current_pos) - 15))))
    if(length(contig) == 1)
      hit_table$contig[row] <- contig
  }
  
  return(hit_table)  
  
}


test_hits_function <- function(blast_csv, ideal_hit_length, contig_dir_loc){
  #browser()
  left_end_isos <- dplyr::count(blast_csv, subject)
  blast_csv <- blast_csv %>% mutate(ori = ifelse(sstart < send, "forward","reverse"))
  hits_df <- NULL
  tot_length <- nrow(left_end_isos)
  total_length <- nchar(as.character(tot_length))
  cat("\n")
  for(row in 1:nrow(left_end_isos)){
    cat(paste("\r On row:", formatC(row, width = total_length, flag = "0"), "of:", as.character(nrow(left_end_isos))))
    current_iso <- left_end_isos[row,1]
    narrowed_blast <- blast_csv[blast_csv$subject == current_iso,]
    narrowed_blast$index <- seq(1, nrow(narrowed_blast))
    match_frame <- NULL
    
    ## Now we need to add in the contig bounds incorporation. 
    narrowed_blast <- contig_location_adder(current_iso, contig_dir_loc, narrowed_blast, row)
    exact_matches <- filter(narrowed_blast, align >= ideal_hit_length - 10)
    if(nrow(exact_matches) > 0){
      nums <- seq(1,nrow(exact_matches))
      extra_names <- paste(exact_matches$subject[1], "_", as.character(nums))
      
      match_frame <- exact_matches %>% mutate(hit_id = extra_names)
      narrowed_blast <- narrowed_blast %>% filter(!(index %in% match_frame$index))
      
    }
    if(nrow(narrowed_blast) > 1){
      ## Now we look to combine hits 
      #browser()
      narrowed_blast <- narrowed_blast %>% rowwise() %>% mutate(sval = min(c(sstart,send))) %>% as.data.frame()
      narrowed_blast <- narrowed_blast[order(narrowed_blast$sval),]
      n <- nrow(narrowed_blast)
      merger <- 2
      start_iso <- NULL
      while(merger <= n){
        if(is.null(nrow(start_iso)))
          start_iso <- narrowed_blast[(merger - 1),]
        next_isos <- narrowed_blast[merger:n,]
        
        if(start_iso$ori == "forward"){
          next_isos <- next_isos %>% filter(ori == "forward") %>% filter(contig == start_iso$contig) %>% 
            filter(send <= start_iso$send + 2000) %>% filter(qstart >= (start_iso$qend - 200)) %>%
            filter(qstart <= (start_iso$qend + 200))
          if(nrow(next_isos) > 0){
            #browser()
            start_iso <- start_iso %>% mutate(qend = next_isos[1,"qend"]) %>% mutate(send = next_isos[1,"send"]) %>%
              mutate(align = qend - qstart) 
            if(start_iso$align >= ideal_hit_length - 10){
              match_frame <- bind_rows(match_frame, start_iso)
              start_iso <- NULL
            }
          }else{
           start_iso <- NULL 
          }
        }else{
          next_isos <- next_isos %>% filter(ori == "reverse") %>% filter(contig == start_iso$contig) %>% 
            filter(send <= start_iso$sstart + 2000) %>% filter(qend >= (start_iso$qstart - 200)) %>%
            filter(qend <= (start_iso$qstart + 200))
          if(nrow(next_isos) > 0){
            start_iso <- start_iso %>% mutate(qstart = next_isos[1,"qstart"]) %>% mutate(sstart = next_isos[1,"sstart"]) %>%
              mutate(align = qend - qstart) 
            if(start_iso$align >= ideal_hit_length - 10){
              match_frame <- bind_rows(match_frame, start_iso)
              start_iso <- NULL
            }
          }else{
            start_iso <- NULL 
          }
        }
        
        merger <- merger + 1
        }
      }
      
      
      
      hits_df <- bind_rows(hits_df, match_frame)
      
    }
  
  cat("\n")
  return(hits_df)    
    
  }

## Main run 

input_args <- get_input()
left_end_csv <- read.csv(input_args[1],
                         header = FALSE, stringsAsFactors = FALSE)
blast_cols <- c("query","subject","pid","align","gap","mismatch","qstart","qend", "sstart","send","eval","bitscore")

colnames(left_end_csv) <- blast_cols
cat("\n","Merging left hits")
test_left <- test_hits_function(left_end_csv,2891, input_args[2])

write.csv(test_left, file = "./left_end_merged.csv", row.names = FALSE, quote = FALSE)
write.csv(test_right, file = "./right_end_merged.csv", row.names = FALSE, quote = FALSE)




