stan_model_cache <- function(name, dir = 'stan') {
  # Filenames
  bfile <- file.path(dir, name)
  code_file <- paste0(bfile, '.stan')
  hash_file <- paste0(bfile, '.hash')
  model_file <- paste0(bfile, '.rds')
  
  # Read source code
  if (!file.exists(code_file)) 
    stop(sprinf('No file named %s', code_file))
  
  # Hash source code
  stan_source <- read_lines(code_file) %>% 
    sub('\\s*(//.*)?$', '', .)
  new_hash <- digest::digest(stan_source)
  
  # What to do if we need to update model
  recompile <- function() {
    message('Compiling stan model from source...')
    model <- stan_model(file = code_file)
    write_lines(new_hash, hash_file)
    saveRDS(model, model_file)
    return(model)
  }
  
  # Read in previous hash to compare 
  if (file.exists(hash_file)) {
    old_hash <- read_lines(hash_file)
  } else {
    message(sprintf("Can't find hash file at %s", hash_file))
    return(recompile())
  }
  
  if (!file.exists(model_file) ||
      new_hash != old_hash) {
    return(recompile())
  } else {
    # otherwise return already compiled model
    message(sprintf('Loading saved stan model %s...', name))
    model <- read_rds(model_file)
    return(model)
  }
}