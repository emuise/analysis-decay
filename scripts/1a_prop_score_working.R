library(terra)
library(tidyverse)
library(tidyterra)
library(data.table)
library(arrow)
library(MatchIt)
library(budR)
source(here::here("scripts", "0a_make_folders.R"))
my_theme <- theme_bw() +
  theme(panel.grid = element_blank())

theme_set(my_theme)

forests <- rast(here::here(raster_loc, "forests.dat"))

bcb <- vect(here::here("data", "shapefiles", "bcb.shp"))

df_save <- function(rast) {
  name <- names(rast)
  layer_sname <- here::here(zone_dir, paste0(name, ".parquet"))
  print(layer_sname)
  
  if (!file.exists(layer_sname)) {
    df <- as.data.frame(rast, xy = T)
    write_parquet(df, sink = layer_sname)
  }
  
  
  return()
}

combine_csvs <- function(folder) {
  save_name <-
    here::here(dirname(folder), paste0(basename(folder), ".parquet"))
  
  print(save_name)
  
  if (!file.exists(save_name)) {
    csvs <- list.files(folder, full.names = T)
    
    joined_csv <- map(csvs, read_parquet) %>%
      map(mutate, y = round(y, digits = 1)) %>%
      reduce(.x = ., .f = dplyr::left_join) %>%
      tibble()
    
    write_parquet(joined_csv, sink = save_name)
  }
}

count_na <- function(parq_loc) {
  if ("character" %in% class(parq_loc)) {
    df <- read_parquet(parq_loc)
    zone <- parq_loc %>%
      tools::file_path_sans_ext() %>%
      basename()
    
    df["zone"] = zone
  } else {
    df <- parq_loc
  }
  zone <- df %>% pull(zone) %>% unique()
  
  print(glue::glue("Count NA: {zone}"))
  
  df %>% summarise_all(~ sum(is.na(.))) %>%
    mutate(zone = zone) %>%
    relocate(zone)
}

# generate structure csvs, don't join them
dat_info <-
  list.files(
    raster_vars,
    pattern = ".dat$",
    full.names = T,
    recursive = T
  ) %>%
  map_dfr(file.info) %>%
  arrange(size) %>%
  rownames_to_column(var = "name") %>%
  tibble()

ordered_rasts <- dat_info %>%
  pull(name)

# new loop through every layer
for (path in ordered_rasts) {
  zone <- path %>%
    dirname() %>%
    basename()
  
  var <- path %>%
    basename() %>%
    tools::file_path_sans_ext()
  
  zone_dir <- here::here(csv_vars, zone)
  
  dir.create(zone_dir, showWarnings = F)
  
  layer_sname <- path %>%
    basename() %>%
    tools::file_path_sans_ext() %>%
    paste0(., ".parquet")
  
  save_name <- here::here(csv_vars, zone, layer_sname)
  
  if (!file.exists(save_name)) {
    print(glue::glue("{zone} {var}"))
    r <- rast(path)
    names(r) <- var
    df <- as.data.frame(r, xy = T)
    write_parquet(df, sink = save_name)
  }
  
}

var_folds <- list.dirs(csv_vars, full.names = T, recursive = F)

map(var_folds, combine_csvs, .progress = T)

var_parqs <-
  list.files(csv_vars, pattern = ".parquet$", full.names = T)

# var_nas <- map_dfr(var_parqs, count_na, .progress = T)

# var_nas

# generate covariate csvs and join them
dat_info <-
  list.files(raster_cov,
             pattern = ".dat$",
             full.names = T) %>%
  map_dfr(file.info) %>%
  arrange(size) %>%
  rownames_to_column(var = "name") %>%
  tibble()

ordered_rasts <- dat_info %>%
  pull(name)

# new loop through every layer
for (path in ordered_rasts) {
  zone <- path %>%
    tools::file_path_sans_ext() %>%
    basename()
  print(zone)
  
  zone_dir <- here::here(csv_cov, zone)
  
  dir.create(zone_dir, showWarnings = F)
  
  r <- rast(path)
  
  rasts <- as.list(r)
  
  map(rasts, df_save, .progress = T)
}

cov_folds <- list.dirs(csv_cov, full.names = T, recursive = F)

map(cov_folds, combine_csvs, .progress = T)

parqs <- list.files(csv_cov, pattern = ".parquet", full.names = T)

# cov_nas <- map_dfr(parqs, count_na, .progress = T)

# cov_nas

load_parq <- function(parq) {
  name <- basename(parq) %>%
    tools::file_path_sans_ext()
  
  parq %>%
    open_dataset() %>%
    mutate(pa = ifelse(is.na(pa), 0, 1)) %>%
    rename(class_val = forests) %>%
    left_join(keys$vlce %>%
                select(class_val, class_name)) %>%
    select(-class_val) %>%
    #drop_na() %>% # 2 pixels in cwh that are going
    mutate(zone = name) %>%
    collect()
}

create_treatment <- function(df) {
  df %>%
    mutate(treat = bc_albers_footprint < 4 * pa) %>%
    relocate(zone, treat, pa, class_name, bc_albers_footprint)
}

create_bins <- function(df, bins = 4) {
  df %>%
    select(-city_distance,-pop_count,-pop_density,-travel_time_to_cities_bc) %>%
    mutate(across(clim_MAP:slope, \(x) ntile(x, n = bins), .names = "{.col}_bin"))
}

tile_names <-
  c("clim_MAP",
    "clim_MAT",
    "clim_MCMT",
    "clim_MWMT",
    "DEM",
    "slope")

x <- parqs[4]
y <- var_parqs[4]

joined_paths <- map2_vec(
  .x = parqs,
  .y = var_parqs,
  .f = \(x, y) {
    zone1 <- x %>% tools::file_path_sans_ext() %>% basename()
    zone2 <- y %>% tools::file_path_sans_ext() %>% basename()
    
    if (zone1 != zone2) {
      break
    }
    
    print(zone1)
    
    save_loc <- here::here(csv_all, paste0(zone1, ".parquet"))
    
    if (!file.exists(save_loc)) {
      parq <- load_parq(x) %>%
        create_treatment() %>%
        create_bins() %>%
        select(-clim_MAP, -clim_MAT, -clim_MCMT, -clim_MWMT, -DEM, -slope) %>%
        unite(col = "tiles", ends_with("bin"), sep = "-") %>%
        as.data.table()
      
      var_parq <- read_parquet(y) %>%
        as.data.table()
      
      print(glue::glue("covariate rows: {nrow(parq)}"))
      print(glue::glue("variate rows: {nrow(var_parq)}"))
      
      joined <- parq[var_parq, on = c("x", "y")]
      
      print(glue::glue(
        "OG Dupes {var_parq %>% select(-x, -y) %>% duplicated() %>% sum()}"
      ))
      #
      print(
        glue::glue(
          "Joined Dupes {joined %>% select(-names(parq)) %>% duplicated() %>% sum()}"
        )
      )
      
      #
      # joined <-
      #   merge(as.data.table(parq),
      #         as.data.table(var_parq),
      #         by = c("x", "y"),
      #         all = T)
      
      print(glue::glue("joined rows: {nrow(joined)}"))
      
      print(count_na(joined))
      
      write_parquet(joined, save_loc)
      
      
    }
    gc()
    save_loc
  },
  .progress = T
)



# split dataframes

map_vec(joined_paths, \(x) {
  df <- read_parquet(x)
  
  zone <- x %>%
    basename() %>%
    tools::file_path_sans_ext()
  print(zone)
  
  zone_loc <- here::here(split_df_loc, zone)
  dir.create(zone_loc, showWarnings = F)
  
  df %>%
    group_by(tiles, class_name) %>%
    group_split() %>%
    map(.f = \(y) {
      tilename <- y %>%
        pull(tiles) %>%
        unique()
      
      class_name <- y %>%
        pull(class_name) %>%
        unique()
      
      save_name <-
        here::here(zone_loc, glue::glue("{class_name}_{tilename}.parquet"))
      #print(save_name)
      if (!file.exists(save_name)) {
        write_parquet(y, sink = save_name)
      }
      save_name
    })
})

split_files <- list.files(split_df_loc, recursive = T, full.names = T)


# bootstrapping confidence intervals

library(boot)

mean_fun <- function(data, indices) {
  col_means <- sapply(data[indices, ], mean)
  return(col_means)
}

x <- split_files[[5]]

map_dfr(split_files, \(x) {
  df <- read_parquet(x)
  zone <- x %>% dirname() %>% basename()
  tile <- x %>% basename() %>% tools::file_path_sans_ext()
  
  treat <- df %>% filter(treat) %>%
    select(basal_area:VarDHI, -MinDHI)
  
  R <- 1000
  
  # Perform bootstrapping
  boot_results <- boot(data = treat, statistic = mean_fun, R = R)
  
  get_ci <- function(col_index) {
    df <- boot.ci(boot_results, type = "bca", index = col_index)$bca |>
      as.data.frame()
    
    names(df) <- c("conf", "bias", "skewness", "lower", "upper")
    
    df["col_index"] = col_index
    
    df
  }
  
  df_list <- lapply(1:ncol(treat), get_ci)
  
  do.call(bind_rows, df_list)
})

# make table for NN

tiles_nn <- map_dfr(joined_paths, \(x) {
  zone <- x %>%
    basename() %>%
    tools::file_path_sans_ext()
  
  df <- open_dataset(x) %>%
    count(class_name, tiles, treat) %>%
    collect() %>%
    drop_na() %>%
    complete(class_name, tiles, treat, fill = list(n = 0)) %>%
    mutate(treat = ifelse(treat, "t", "ut")) %>%
    pivot_wider(names_from = treat, values_from = n) %>%
    arrange(t) %>%
    separate(tiles, into = tile_names) %>%
    mutate(across(clim_MAP:slope, as.numeric),
           zone = zone) %>%
    relocate(zone)
}, .progress = T)

match_treats <- function(u_split, t_split) {
  nn <- FNN::get.knnx(t_split %>%
                        select(clim_MAP:slope),
                      u_split %>%
                        select(clim_MAP:slope),
                      k = 1)
  
  t_split <- t_split %>%
    mutate(row_no = row_number()) %>%
    select(-class_name, -n)
  
  u_split %>%
    mutate(nn_t = as.numeric(nn$nn.index),
           nn_dist = as.numeric(nn$nn.dist)) %>%
    left_join(t_split,
              by = c("nn_t" = "row_no"),
              suffix = c(".u", ".t"))
  
  
}

tiles_nn %>%
  group_by(zone) %>%
  group_split() %>%
  map(., .f = \(df) {
    t_splits <- df %>%
      filter(t > 3)
    
    ut_splits <- df %>%
      filter(ut != 0, t <= 3)
    
    ut_splits %>%
      group_split(class_name) %>%
      map(\(y) {
        class <- y %>%
          pull(class_name) %>%
          unique()
        
        t <- t_splits %>%
          filter(class_name == class)
        
        if(nrow(t) == 0) {
          t <- t_splits
        }
        
        nn <- FNN::get.knnx(t %>% 
                              select(clim_MAP:slope), 
                            y %>% 
                              select(clim_MAP:slope), 
                            k = 1)
        
        t <- t %>%
          mutate(row_no = row_number())
        
        y_matched <- y %>%
          mutate(nn_ind = as.numeric(nn$nn.index),
                 nn_dist = as.numeric(nn$nn.dist))%>%
          left_join(t,
                    by = c("nn_ind" = "row_no"),
                    suffix = c(".u", ".t"))
      })
    
  })



df <- parqs[[4]] %>%
  load_parq() %>%
  create_treatment() %>%
  create_bins()

unmatched_splits <- df %>%
  filter(t == 0, ut != 0) %>%
  separate(tiles,
           into = c("MAP", "MAT", "MCMT", "MWMT", "elev", "slope"),
           sep = "-") %>%
  mutate(across(MAP:slope, as.numeric)) %>%
  group_by(class_name, zone) %>%
  group_split()

treated_splits <- df %>%
  filter(t > 0) %>%
  separate(tiles,
           into = c("MAP", "MAT", "MCMT", "MWMT", "elev", "slope"),
           sep = "-") %>%
  mutate(across(MAP:slope, as.numeric)) %>%
  group_by(class_name) %>%
  group_split()




treated <- zone_match %>%
  filter(pa == 1, n != 0) %>%
  separate(tiles,
           into = c("MAP", "MAT", "MCMT", "MWMT", "elev", "slope"),
           sep = "-")

treated_splits <- treated %>%
  mutate(across(MAP:slope, as.numeric)) %>%
  group_by(class_name) %>%
  group_split()



nn_match <-
  map2_dfr(.x = unmatched_splits, .y = treated_splits, .f = match_treats)

# plot of matched vs unmatched as a function of nn distance



zone_counts <- zone_match %>%
  pivot_wider(names_from = pa,
              values_from = n) %>%
  rename(ua = 3, pa = 4)

matched_ua <- zone_counts %>%
  filter(ua != 0,
         pa != 0) %>%
  select(class_name, ua) %>%
  mutate(nn_dist = 0)

to_bind_ua <- nn_match %>%
  select(class_name, ua = ut, nn_dist)


bind_rows(matched_ua, to_bind_ua) %>%
  group_by(class_name, nn_dist) %>%
  summarize(n = sum(ua)) %>%
  mutate(per = n / sum(n)) %>%
  ggplot(aes(x = nn_dist, y = per)) +
  geom_line(aes(col = class_name))
# we assess the distance decay of structural and functional metrics in forest
# stands across the forested ecosystems of british columbia

# i want to figure out the amount of unmatched pixels
# to see if it's a large amount or not



matched <- zone_match %>%
  pivot_wider(names_from = pa, values_from = n) %>%
  rename(pa = `1`, ua = `0`) %>%
  filter(pa != 0 & ua != 0)

zone_match %>%
  filter(tiles %in% unmatched & n != 0)

cem_out <-
  matchit(
    treat ~ clim_MAP + clim_MAT + clim_MCMT + clim_MWMT + DEM + slope,
    method = "cem",
    data = zone_treat
  )

cem_out <-
  matchit(
    pa ~ bc_albers_footprint + city_distance + clim_MAP + clim_MAT + clim_MCMT + clim_MWMT + DEM + slope + pop_count + pop_density + travel_time_to_cities_bc,
    data = zone,
    cutpoints = list(
      bc_albers_footprint = bins,
      city_distance = bins,
      clim_MAP = bins,
      clim_MAT = bins,
      clim_MCMT = bins,
      clim_MWMT = bins,
      DEM = bins,
      slope = bins,
      pop_count = bins,
      pop_density = bins,
      travel_time_to_cities_bc = bins
    ),
    method = "cem"
  )

my_match <- zone %>%
  #slice_sample(n = 10000) %>%
  relocate(pa) %>%
  mutate(across(
    bc_albers_footprint:travel_time_to_cities_bc,
    \(x) ntile(x, bins)
  )) %>%
  group_by(
    bc_albers_footprint,
    city_distance,
    clim_MAP,
    clim_MAT,
    clim_MCMT,
    clim_MWMT,
    DEM,
    pop_count,
    pop_density,
    slope,
    travel_time_to_cities_bc,
    class_name
  ) %>% count(pa) %>%
  ungroup() %>%
  mutate(
    columns = paste(
      bc_albers_footprint,
      city_distance,
      clim_MAP,
      clim_MAT,
      clim_MCMT,
      clim_MWMT,
      DEM,
      pop_count,
      pop_density,
      slope,
      travel_time_to_cities_bc,
      class_name,
      sep = "-"
    )
  ) %>%
  select(columns, pa, n) %>%
  complete(columns, pa, fill = list(n = 0))

unmatched = my_match %>%
  filter(n == 0) %>%
  pull(columns)

my_match %>%
  filter(columns %in% unmatched)

my_match %>% count(n == 0)

cem_out <-
  matchit(
    pa ~ bc_albers_footprint + class_name + city_distance + clim_MAP + clim_MAT + clim_MCMT + clim_MWMT + DEM + slope + pop_count + pop_density + travel_time_to_cities_bc,
    data = zone %>%
      slice_sample(n = 10000),
    cutpoints = list(
      bc_albers_footprint = bins,
      city_distance = bins,
      clim_MAP = bins,
      clim_MAT = bins,
      clim_MCMT = bins,
      clim_MWMT = bins,
      DEM = bins,
      slope = bins,
      pop_count = bins,
      pop_density = bins,
      travel_time_to_cities_bc = bins
    ),
    method = "cem"
  )


zone %>%
  slice_sample(n = 100000) %>%
  relocate(pa) %>%
  mutate(pa = ifelse(pa == 1, "PA", "UA")) %>%
  pivot_longer(cols = bc_albers_footprint:travel_time_to_cities_bc) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(col = pa)) +
  facet_wrap(~ name,
             scales = "free")

cem_out <-
  matchit(
    pa ~ bc_albers_footprint + city_distance + clim_MAP + clim_MAT + clim_MCMT + clim_MWMT + DEM + slope + pop_count + pop_density + travel_time_to_cities_bc,
    data = zone,
    cutpoints = list(
      bc_albers_footprint = bins,
      city_distance = bins,
      clim_MAP = bins,
      clim_MAT = bins,
      clim_MCMT = bins,
      clim_MWMT = bins,
      DEM = bins,
      slope = bins,
      pop_count = bins,
      pop_density = bins,
      travel_time_to_cities_bc = bins
    ),
    method = "cem"
  )


#####


quantile_df <- function(x, probs = seq(0, 1, .01)) {
  tibble(val = quantile(x, probs, na.rm = TRUE),
         quant = probs)
}


quants <- zone %>%
  relocate(pa) %>%
  reframe(across(bc_albers_footprint:slope, quantile_df, .unpack = T),
          .by = c(pa, forests)) %>%
  mutate(quantile = bc_albers_footprint_quant) %>%
  select(!ends_with("quant"))

quants %>%
  select(quantile, pa, bc_albers_footprint_val) %>%
  pivot_wider(names_from = pa, values_from = bc_albers_footprint_val) %>%
  ggplot(aes(x = PA, y = UA)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              lty = "dotted")

quants %>%
  pivot_longer(cols = c(-pa, -quantile, -forests)) %>%
  pivot_wider(names_from = pa, values_from = value) %>%
  ggplot(aes(x = PA, y = UA, col = forests), alpha = 0.5) +
  geom_point() +
  geom_abline(
    slope = 1,
    intercept = 0,
    lty = "dotted",
    colour = "red"
  ) +
  facet_wrap(~ name, scales = "free")
