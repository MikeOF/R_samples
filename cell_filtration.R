# Author: Michael Finlayson
#
# Description:
# This class encapsulates the process of filtering cells in a scRNA-seq sample.
# The goal here is to make standardized cell filtration decisions automatic.
# Parameters can be set to enforce hard filtration thresholds when automatic
# thresholds appear to be insufficient.
#

######## -----------------------------------------------------------------------
########                 Parameters
######## ---------------------------------------------------------------

CFParameters = setClass('CFParameters',
  slots = c(
   min_nUMI = 'numeric',
   max_nUMI = 'numeric',
   min_nGene = 'numeric',
   max_nGene = 'numeric',
   min_mt_prop = 'numeric',
   max_mt_prop = 'numeric'
  )
)

setMethod('initialize', 'CFParameters',
function(.Object,
         min_nUMI = -Inf, max_nUMI = Inf,
         min_nGene = -Inf, max_nGene = Inf,
         min_mt_prop = -Inf, max_mt_prop = Inf) {
  .Object@min_nUMI = min_nUMI
  .Object@max_nUMI = max_nUMI
  .Object@min_nGene = min_nGene
  .Object@max_nGene = max_nGene
  .Object@min_mt_prop = min_mt_prop
  .Object@max_mt_prop = max_mt_prop
  return(.Object)
})

######## -----------------------------------------------------------------------
########                 Definition
######## ---------------------------------------------------------------

CellFiltration = setClass('CellFiltration',
  slots = c(
    parameters = 'CFParameters',
    filtration_by_sample = 'list',
    cells_to_keep = 'character'
  )
)

######## -----------------------------------------------------------------------
########                 Initialization
######## ---------------------------------------------------------------

setMethod('initialize', 'CellFiltration',
function(.Object, scdata, parameters = CFParameters()) {
  return(.cf_init(.Object, scdata, parameters = parameters))
})

getCellFiltration = function(scdata, parameters = CFParameters()) {

  # check input types
  stopifnot(is(scdata, 'SCData'))
  stopifnot(is(parameters, 'CFParameters'))

  return(CellFiltration(scdata, parameters))
}

######## -----------------------------------------------------------------------
########                 Private Functions
######## ---------------------------------------------------------------

###### -------------------------------------
######     .cf_init

.cf_init = function(cf, scdata, parameters) {

  # check input types
  stopifnot(is(cf, 'CellFiltration'))
  stopifnot(is(scdata, 'SCData'))
  stopifnot(is(parameters, 'CFParameters'))

  # check inputs
  stopifnot(all(c('sample', 'nUMI', 'nGene', 'mt_prop') %in% colnames(scdata@metadata)))

  metadata = scdata@metadata

  # create a "filtration" for each sample, ie filter each sample
  filtration_by_sample = list()
  for (sample_name in unique(metadata$sample)) {

    # initialize the filtration data with quality metrics
    filtration_data = metadata[metadata$sample == sample_name, c('nGene', 'nUMI', 'mt_prop'), drop=F]

    # filter based on minimums and maximums in the parameters
    filtration_data = .cf_apply_parameter_limits(filtration_data, parameters)

    # first thing, mt_prop
    filtration_data = .cf_get_mt_prop_probablities(filtration_data)
    filtration_data = .cf_apply_mt_prop_prob_limits(filtration_data)

    # Now, nUMI
    filtration_data = .cf_get_nUMI_probablities(filtration_data)
    filtration_data = .cf_apply_nUMI_prob_limits(filtration_data)

    # Now, nGene
    filtration_data = .cf_get_nGene_probabilities(filtration_data)
    filtration_data = .cf_apply_nGene_prob_limits(filtration_data)

    # get cells to keep by applying the cutoffs
    filtration_data = .cf_get_cells_to_keep(filtration_data)

    filtration_by_sample[[sample_name]] = filtration_data
  }

  # compose cell filtration object
  cf@filtration_by_sample = filtration_by_sample

  for(sample_name in names(filtration_by_sample)) {
    sample_filtration_data = filtration_by_sample[[sample_name]]
    sample_cells_to_keep = rownames(sample_filtration_data)[sample_filtration_data$within_all_limits]
    cf@cells_to_keep = union(cf@cells_to_keep, sample_cells_to_keep)
  }

  cf@parameters = parameters

  return(cf)
}

###### -------------------------------------
######     .cf_apply_parameter_limits

.cf_apply_parameter_limits = function(filtration_data, parameters) {

  # determine which cells are outside thresholds
  gt_max_nUMI = filtration_data$nUMI > parameters@max_nUMI
  lt_min_nUMI = filtration_data$nUMI < parameters@min_nUMI
  gt_max_nGene = filtration_data$nGene > parameters@max_nGene
  lt_min_nGene = filtration_data$nGene < parameters@min_nGene
  gt_max_mt_prop = filtration_data$mt_prop > parameters@max_mt_prop
  lt_min_mt_prop = filtration_data$mt_prop < parameters@min_mt_prop

  # determine which cells are filtered "beyond_limits"
  beyond_limits = gt_max_nUMI | lt_min_nUMI | gt_max_nGene
  beyond_limits = beyond_limits | lt_min_nGene | gt_max_mt_prop
  beyond_limits = beyond_limits | lt_min_mt_prop

  # create the filtration data
  filtration_data = data.frame(filtration_data,
                               gt_max_nUMI = gt_max_nUMI,
                               lt_min_nUMI = lt_min_nUMI,
                               gt_max_nGene = gt_max_nGene,
                               lt_min_nGene = lt_min_nGene,
                               gt_max_mt_prop = gt_max_mt_prop,
                               lt_min_mt_prop = lt_min_mt_prop,
                               within_parameter_limits = !beyond_limits)

  return(filtration_data)
}

###### -------------------------------------
######     .cf_get_mt_prop_probablities

.cf_get_mt_prop_probablities = function(filtration_data) {

  # determine which cells to fit on by taking higher quality cells
  # according to median nUMI and nGene values
  log_nUMI = log1p(filtration_data$nUMI)
  log_nUMI_median = median(log_nUMI)
  log_nUMI_sd = sqrt( sum( (log_nUMI - log_nUMI_median)^2 ) / (length(log_nUMI) - 1))

  nGene = filtration_data$nGene
  nGene_fit = fitdistrplus::fitdist(nGene, distr = 'norm')
  nGene_mean = nGene_fit$estimate['mean']
  nGene_sd = nGene_fit$estimate['sd']

  fit_data_select = abs(log_nUMI - log_nUMI_median) < log_nUMI_sd & abs(nGene - nGene_mean) < nGene_sd

  # log transform the data
  log_mt_prop = log1p(filtration_data$mt_prop)
  names(log_mt_prop) = rownames(filtration_data)

  # get data to fit on
  fit_data = log_mt_prop[fit_data_select]

  # fit the data
  fit_mean = median(fit_data)
  fit_sd = sqrt( sum( (fit_data - fit_mean)^2 ) / length(fit_data) )

  # compute the probabilities
  mt_prop_prob = pnorm(log_mt_prop, mean = fit_mean, sd = fit_sd)

  filtration_data = data.frame(filtration_data, mt_prop_prob=mt_prop_prob)

  return(filtration_data)
}

###### -------------------------------------
######     .cf_apply_mt_prop_prob_limits

.cf_apply_mt_prop_prob_limits = function(filtration_data) {

  gt_mt_prop_max_prob_limit = filtration_data$mt_prop_prob > 0.96
  within_mt_prop_prob_limits = !(gt_mt_prop_max_prob_limit)

  filtration_data = data.frame(filtration_data,
                               gt_mt_prop_max_prob_limit=gt_mt_prop_max_prob_limit,
                               within_mt_prop_prob_limits=within_mt_prop_prob_limits)

  return(filtration_data)
}

###### -----------------------------------
######     .cf_get_nUMI_probablities

.cf_get_nUMI_probablities = function(filtration_data) {

  # determine which cells to fit on
  fit_data_select = filtration_data$within_mt_prop_prob_limits

  # log transform the data
  log_nUMI = log1p(filtration_data$nUMI)

  # get data to fit on
  fit_data = log_nUMI[fit_data_select]

  # fit the data
  median_log_nUMI = median(fit_data)
  sd_log_nUMI = sqrt( sum( (fit_data - median_log_nUMI)^2) / (length(log_nUMI) - 1) )

  # compute the probabilities
  nUMI_prob = pnorm(log_nUMI, mean = median_log_nUMI, sd = sd_log_nUMI)

  filtration_data = data.frame(filtration_data, nUMI_prob=nUMI_prob)
  return(filtration_data)
}

###### -----------------------------------
######     .cf_apply_nUMI_prob_limits

.cf_apply_nUMI_prob_limits = function(filtration_data) {

  gt_nUMI_max_prob_limit = filtration_data$nUMI_prob > 0.98
  lt_nUMI_min_prob_limit = filtration_data$nUMI_prob < 0.02
  within_nUMI_prob_limits = !(gt_nUMI_max_prob_limit | lt_nUMI_min_prob_limit)

  filtration_data = data.frame(filtration_data,
                               gt_nUMI_max_prob_limit=gt_nUMI_max_prob_limit,
                               lt_nUMI_min_prob_limit=lt_nUMI_min_prob_limit,
                               within_nUMI_prob_limits=within_nUMI_prob_limits)

  return(filtration_data)
}

###### -----------------------------------
######     .cf_get_nGene_probabilities

.cf_get_nGene_probabilities = function(filtration_data) {

  # determine which cells to fit on
  fit_data_select = filtration_data$within_mt_prop_prob_limits

  nGene = filtration_data$nGene

  # get data to fit on
  fit_data = nGene[fit_data_select]

  # fit the data
  dist_fit = MASS::fitdistr(fit_data, densfun = 'normal')

  # compute probabilities
  nGene_prob = pnorm(nGene,
                     mean = dist_fit$estimate['mean'],
                     sd = dist_fit$estimate['sd'])

  filtration_data = data.frame(filtration_data, nGene_prob=nGene_prob)
  return(filtration_data)
}

###### -----------------------------------
######     .cf_apply_nGene_prob_limits

.cf_apply_nGene_prob_limits = function(filtration_data) {

  gt_nGene_max_prob_limit = filtration_data$nUMI_prob > 0.99
  lt_nGene_min_prob_limit = filtration_data$nUMI_prob < 0.01
  within_nGene_prob_limits = !(gt_nGene_max_prob_limit | lt_nGene_min_prob_limit)

  filtration_data = data.frame(filtration_data,
                               gt_nGene_max_prob_limit=gt_nGene_max_prob_limit,
                               lt_nGene_min_prob_limit=lt_nGene_min_prob_limit,
                               within_nGene_prob_limits=within_nGene_prob_limits)

  return(filtration_data)
}

###### -----------------------------
######     .cf_get_cells_to_keep

.cf_get_cells_to_keep = function(filtration_data) {

  within_all_limits = filtration_data$within_parameter_limits
  within_all_limits = filtration_data$within_mt_prop_prob_limits & within_all_limits
  within_all_limits = filtration_data$within_nUMI_prob_limits & within_all_limits
  within_all_limits = filtration_data$within_nGene_prob_limits & within_all_limits

  filtration_data = data.frame(filtration_data, within_all_limits=within_all_limits)
  return(filtration_data)
}

######## -----------------------------------------------------------------------
########                 Public Functions
######## ---------------------------------------------------------------

###### --------------------
######     cf_nullify

cf_nullify = function(cf) {

  # check input type
  stopifnot(is(cf, 'CellFiltration'))

  cf@parameters = CFParameters()
  cf@filtration_by_sample = list()
  cf@cells_to_keep = NA_character_

  return(cf)
}

###### --------------------
######     cf_get_plots

cf_get_plots = function(cf) {

  # check input type
  stopifnot(is(cf, 'CellFiltration'))

  # create a plot for each sample
  plot_list = list()
  for (sample_name in names(cf@filtration_by_sample)) {

    filtration_data = cf@filtration_by_sample[[sample_name]]

    # collect quality metrics for each sample
    qual_metrics = filtration_data[,c('nGene','nUMI','mt_prop'),drop=F]

    # get metrics before filtration
    qual_metrics_before = qual_metrics
    qual_metrics_before = data.frame(time=rep('before', times=nrow(qual_metrics_before)),
                                     qual_metrics_before)

    # get metrics after filtration
    qual_metrics_after = qual_metrics[filtration_data$within_all_limits,,drop=F]
    qual_metrics_after = data.frame(time=rep('after', times=nrow(qual_metrics_after)),
                                    qual_metrics_after)

    # create data to plot with
    pdata = rbind(qual_metrics_before, qual_metrics_after)
    pdata[,'mt_prop'] = pdata[,'mt_prop'] + 0.00001
    pdata[,'nUMI'] = pdata[,'nUMI'] + 0.001

    # determine aesthetics to use
    stroke_size = 1
    alpha_val = .8
    if (nrow(qual_metrics) > 50) {
      stroke_size = .5
      alpha_val = .5
    }
    if (nrow(qual_metrics) > 200) {
      stroke_size = .2
      alpha_val = .2
    }
    if (nrow(qual_metrics) > 1000) {
      stroke_size = .05
      alpha_val = .01
    }
    if (nrow(qual_metrics) > 2000) {
      stroke_size = .01
      alpha_val = .001
    }

    # create a before and after plot for each sample
    plots = list()
    for (metric in colnames(qual_metrics)) {

      ggp = ggplot2::ggplot(
        reshape2::melt(pdata[,c('time',metric)], id.vars = 'time'),
        ggplot2::aes(x=time, y=value)
      )
      ggp = ggp + ggplot2::theme_minimal()
      ggp = ggp + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      )
      ggp = ggp + ggplot2::geom_violin(ggplot2::aes(fill='red'))
      ggp = ggp + ggplot2::geom_jitter(ggplot2::aes(alpha = .2, stroke = stroke_size))
      ggp = ggp + ggplot2::ggtitle(metric)
      ggp = ggp + ggplot2::scale_fill_discrete(guide=F)
      ggp = ggp + ggplot2::scale_alpha(guide=F)

      if (metric %in% c('nUMI', 'mt_prop')) { ggp = ggp + ggplot2::scale_y_continuous(trans='log10') }

      plots[[metric]] = ggp

    }
    filtration_plot = cowplot::plot_grid(plotlist=plots, ncol=3)

    # paste plots for the three metrics together and draw a title
    title = cowplot::ggdraw() + cowplot::draw_label(sample_name, fontface='bold')
    filtration_plot = cowplot::plot_grid(title, filtration_plot, ncol=1, rel_heights=c(0.1, 1))

    plot_list[[sample_name]] = filtration_plot
  }

  return(plot_list)
}

###### ---------------
######     cf_save

cf_save = function(cf, save_dir, dev_to_use = 'pdf') {

  # check input types
  stopifnot(is(cf, 'CellFiltration'))
  stopifnot(is(save_dir, 'character'))
  stopifnot(is(dev_to_use, 'character'))

  # work out that trailing slash ambiguity
  save_dir = paste(strsplit(save_dir, '/')[[1]], collapse='/')

  # don't overwrite with files or directories
  stopifnot(!dir.exists(save_dir) && !(file.exists(save_dir)))
  dir.create(save_dir)

  # save plots and filtration data for each sample
  plot_by_sample_list = cf_get_plots(cf)
  for (sample in names(plot_by_sample_list)) {

    sample_dir = paste(save_dir, sample, sep = '/')
    dir.create(sample_dir)

    # save the filtration data
    fdata_filename = paste(sample, 'filtration_data.csv', sep = '_')
    fdata_filename = paste(sample_dir, fdata_filename, sep = '/')

    filtration_data = cf@filtration_by_sample[[sample]]
    filtration_data = data.frame(cell_name=rownames(filtration_data), filtration_data)

    write.csv(x = filtration_data, file = fdata_filename, row.names = F)

    # save the plot
    plot_filename = paste(sample, paste('_filtration', dev_to_use, sep = '.'), sep = '')
    plot_filename = paste(sample_dir, plot_filename, sep = '/')

    cf_plot = plot_by_sample_list[[sample]]

    ggplot2::ggsave(filename = plot_filename, plot = cf_plot, device = dev_to_use, height = 7, width = 9)
  }
}
