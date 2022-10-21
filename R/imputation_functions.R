
#' Preparation of input files for imputation
#'
#' This function is used as a nested function in phcp_impute()
#' @param haps Phased reference genomes in .haps format, which is the output of SHAPEIT2.
#' @param genetic_map genetic map for the reference genome with 3 columns: chromosome number, genetic distance in Morgans, physical position. The columns header should be: chr, genetic_distance, bp
#' @param ancient_tped pseudo-haplotype to impute in Plink tped format
#' @param chr the number of chromosome to impute
#' @return an input for the imputation
#' @seealso [PHCPImpute::phcp_impute] which is the imputation function that uses this function
#' @export
#' @examples
input = function(haps, genetic_map, ancient_tped, chr) {
  colnames(haps)[1:5] = c("chr","id","bp","a0","a1")
  colnames(ancient_tped) = c("chr","id","genetic_distance","bp","ph1","ph2")
  all_genomes <- inner_join(haps, ancient_tped, by = c("chr", "bp")) %>%
    filter( ph1 == a0 | ph1 == a1 | ph1 == 0) %>%
    select(!c("id.x","id.y","genetic_distance","ph2")) %>%
    mutate(ancient_binary = ifelse(ph1 == 0, NA, ifelse(ph1 == a0, 0, 1)), .after = 'a1') %>%
    select(-ph1) %>%
    inner_join(genetic_map,.) %>%
    arrange(chr, bp)
  return(all_genomes)
}

#' Choose most informative donors for imputation
#'
#' The donors are chosen based on K. Wasik et al, 2021
#' This function is used as a nested function in phcp_impute()
#' @param data the input file that was generated with [PHCPImpute::input]
#' @param n_donors the number of donors' haplotypes to choose
#' @param n_info the number of columns in the data file before the columns of the haplotypes
#' @param freq the maximal minor allele frequency in variants to consider when choosing the donors
#' @return the data file with the chosen haplotypes
#' @seealso [PHCPImpute::phcp_impute] which is the imputation function that uses this function
#' @seealso [PHCPImpute::input] that generates the data input for this function
#' @export
#' @examples
#' @references K. Wasik et al., Comparing low-pass sequencing and genotyping for trait mapping in pharmacogenetics. BMC genomics 22, 197 (2021).
filter_donors = function(data, n_donors, n_info, freq){
  data = data %>% mutate(ratio = rowMeans(.[ ,(all_of(n_info) + 1):ncol(.)]), .after = all_of(n_info))
  n_info = n_info + 1
  fT1 = (data$ratio > freq) & (data == 0)
  fT2 = (data$ratio < 1-freq) & (data == 1)
  fT = fT1 + fT2
  
  fM1 = (data$ratio > freq) & (data == 0) & (!is.na(data$ancient_binary) &  data$ancient_binary== 0)
  fM2 = (data$ratio < 1-freq) & (data == 1) & (!is.na(data$ancient_binary) & data$ancient_binary == 1)
  fM = fM1 + fM2
  
  fT = colSums(fT)
  fM = colSums(fM)
  ranks = fM/fT
  ranks[1:n_info] = 1
  ranks[is.nan(ranks)] = 0
  top_donors = order(ranks, decreasing=TRUE)[1:(n_info+n_donors)]
  data <- data[, top_donors]
  output = list(data = data, donors = top_donors[-(1:n_info)]-n_info)
  return(output)
}

#' Imputation for pseudo-haplotypes
#'
#' This function perform imputation of pseudo-haplotypes. The imputation is based on Lee and Li & Stephens hidden Markov model, as implemented in ChromoPainter, but with the hidden state representing a pair of haplotypes from the reference panel.
#' @param haps a data frame of Phased reference genomes in .haps format, which is the output of SHAPEIT2.
#' @param ancient_tped a data frame of pseudo-haplotype to impute in Plink tped format.
#' @param genetic_map a data frame of genetic map for the reference genome with 3 columns: chromosome number, genetic distance in Morgans, physical position. The columns header should be: chr, genetic_distance, bp
#' @param chr the number of chromosome to impute (integer).
#' @param choose_donors the number of haplotypes to choose from the reference data for the imputation. If FALSE, all haplotypes in the data will be used.
#' @param freq_donors the maximal minor allele frequency in variants to consider when choosing the donors
#' @param variants a data frame with variants to impute . If the arguments is missing, all variants of the reference panel will be be imputed.
#' @return A data frames with the posterior probabilities of the possible genotypes in the imputed variants.
#' @export
#' @example
#' @references N. Li, M. Stephens, Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data. Genetics 165, 2213-2233 (2003)
#' @references D. J. Lawson, G. Hellenthal, S. Myers, D. Falush, Inference of population structure using dense haplotype data. PLoS Genet 8, e1002453 (2012)
phcp_impute = function(haps, ancient_tped, genetic_map, chr, choose_donors = F, freq_donors = 0.5, variants) {
  colnames(haps)[1:5] = c("chr","id","bp","a0","a1")
  all_genomes = input(haps, genetic_map, ancient_tped, chr)
  n_info = all_genomes %>% select(!starts_with("V")) %>% ncol(.)
  if(!choose_donors){
    data = all_genomes
  } else{
    filtered_input = filter_donors(data = all_genomes, n_donors = choose_donors, n_info, freq_donors)
    data = filtered_input[["data"]]
  }
  theta <- 0.00138488 # The mutation rate
  N <- 64.5698 # The effective population size
  
  reference_genomes <- as.matrix(data %>% select(starts_with("V")))
  mode(reference_genomes) <- "logical"
  
  
  J <- ncol(reference_genomes) #number of donors haplotypes
  prior <- rep(1/J,J) #vector of a-priori copying probabilities from one haplotype
  prior.prior <- outer(prior,prior)
  ancient <- as.matrix(data$ancient_binary)
  informative_snps <- !is.na(ancient)
  
  #remove variants that are missing in the ancient genome from both donors data and ancient genome
  df <- reference_genomes[informative_snps,]
  ancient <- ancient[informative_snps,]
  
  L <- nrow(df) #number of variants
  
  
  
  ## calculate the gentic distance between each pair of following variants
  rho <- N*c(diff(data$genetic_distance[informative_snps]), Inf)
  
  #forward algorithm (with scaling)
  time = Sys.time()
  
  Theta <- matrix(theta,J,J)
  
  df.snp <- df[1,1:J]
  
  ancient.equal.df <- (ancient[1] == df.snp)*(1-2*theta)
  A.ne.B <- outer(df.snp, df.snp, function(x,y) x != y)
  emission_l <- sweep(Theta, 1, ancient.equal.df, "+")
  emission_l[A.ne.B] <- 0.5
  alpha_hat <- emission_l*prior.prior
  ct <- numeric(L) #scaling vector
  ct[1] <- 1/sum(alpha_hat)
  alpha_hat <- alpha_hat*ct[1]
  prob_array <- array(dim = c(J,J,L))
  prob_array[ , ,1] <- alpha_hat
  for(l in 2:L) {
    expr <- exp(-rho[l-1])
    recomb2 <- (1-expr)^2*(1/(J^2))
    recomb1 <- expr*(1-expr)*(1/J)*rowSums(alpha_hat)
    recomb02 <- exp(-2*rho[l-1])*alpha_hat + recomb2
    recomb021 <- sweep(recomb02, 1, recomb1, "+")
    tot_rec <- sweep(recomb021, 2, recomb1, "+")
    
    df.snp <- df[l,1:J]
    ancient.equal.df <- (ancient[l] == df.snp)*(1-2*theta)
    A.ne.B <- outer(df.snp, df.snp, function(x,y) x != y)
    emission_l <- sweep(Theta, 1, ancient.equal.df, "+")
    emission_l[A.ne.B] <- 0.5
    
    alpha_hat <- emission_l * tot_rec
    ct[l] <- 1/sum(alpha_hat)
    alpha_hat <- ct[l]*alpha_hat
    prob_array[ , ,l] <- alpha_hat
  }
  
  print("finish forward")
  print(Sys.time() - time)
  
  mutations_prob = array(dim = c(J,J,L))
  beta_hat <- ct[L]*matrix(1, nrow = J, ncol = J)
  genotype_prob = matrix(nrow = L, ncol  = 3)
  colnames(genotype_prob) = c("0/0","0/1","1/1")
  states_probs = beta_hat * prob_array[ , ,L] / sum(beta_hat * prob_array[ , ,L])
  alleles <- outer(df[L,1:J], df[L,1:J], "+")
  genotype_prob[L, ] = c(sum(states_probs[alleles == 0]), sum(states_probs[alleles == 1]), sum(states_probs[alleles == 2]))
  mutations_prob[ , ,L] = states_probs
  for(l in (L-1):1) {
    expr <- exp(-rho[l])
    
    df.snp <- df[l+1,1:J]
    ancient.equal.df <- (ancient[l+1] == df.snp)*(1-2*theta)
    A.ne.B <- outer(df.snp, df.snp, function(x,y) x != y)
    emission_l <- sweep(Theta, 1, ancient.equal.df, "+")
    emission_l[A.ne.B] <- 0.5
    
    b <- beta_hat*emission_l
    recomb2 <- (1-expr)^2*(1/(J^2))*sum(b)
    recomb1 <- expr*(1-expr)*(1/J)*rowSums(b)
    recomb02 <- exp(-2*rho[l])*b + recomb2
    recomb021 <- sweep(recomb02, 1, recomb1, "+")
    beta_hat <- sweep(recomb021, 2, recomb1, "+")
    
    beta_hat <- beta_hat*ct[l]
    #####new lines####
    states_probs = beta_hat * prob_array[ , ,l] / sum(beta_hat * prob_array[ , ,l])
    alleles <- outer(df[l,1:J], df[l,1:J], "+")
    genotype_prob[l, ] = c(sum(states_probs[alleles == 0]), sum(states_probs[alleles == 1]), sum(states_probs[alleles == 2]))
    mutations_prob[ , ,l] = states_probs
  }
  print("finished backward")
  
  mutations_info = haps %>%
    inner_join(., genetic_map) %>%
    arrange(bp)
  
  if(!missing("variants")){
    mutations_info = mutations_info %>%
      inner_join(y = variants, by = c("chr", "bp")) %>%
      filter((Ref == a0 & Alt == a1) | (Ref == a1 & Alt == a0)) %>%
      arrange(bp)
  }
  
  mutations_genotype = mutations_info %>% select(starts_with("V"))
  if (choose_donors){
    mutations_genotype = mutations_genotype %>% select(filtered_input[["donors"]])
  }
  mutations_genotype = as.matrix(mutations_genotype)
  mutations_output = matrix(nrow = nrow(mutations_info), ncol = 3)
  colnames(mutations_output) = c("0/0","0/1","1/1")
  for (i in 1:nrow(mutations_info)){
    alleles <- outer(mutations_genotype[i, ], mutations_genotype[i, ], "+")
    closest_snps = which.min(abs(mutations_info$genetic_distance[i] - data$genetic_distance[informative_snps]))
    probs = mutations_prob[,,closest_snps]
    mutations_output[i, ] = c(sum(probs[alleles == 0]), sum(probs[alleles == 1]), sum(probs[alleles == 2]))
  }
  mutations_output = data.frame(as.data.frame(mutations_info)[,c("chr","bp","a0","a1")], mutations_output)
  
  
  
  print("finished imputation")
  print(Sys.time() - start_time)
  
  return(mutations_output)
}



