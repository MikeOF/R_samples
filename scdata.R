# Author: Michael Finlayson
#
# Description:
# This object will create a list of umi_counts and metadata for each library
# that can be found within or below the current directory or directory
# provided
#

######## -----------------------------------------------------------------------
########                 Definition
######## ---------------------------------------------------------------

SCData = setClass('SCData',
  slots=c(
    umi_counts = 'Matrix',
    metadata = 'data.frame'
  )
)

######## -----------------------------------------------------------------------
########                 Initialization
######## ---------------------------------------------------------------

setMethod('initialize', 'SCData', function(.Object, umi_counts, metadata=NULL) {

  if (!is.null(metadata)) { stopifnot(all(colnames(umi_counts)==rownames(metadata))) }

  metadata = .scd_set_count_derived_metadata(umi_counts, metadata)

  .Object@umi_counts = umi_counts
  .Object@metadata = metadata

  return(.Object)
})

getSCData = function(umi_counts, metadata=NULL) {

  # check input types
  stopifnot(is(umi_counts, 'Matrix'))
  stopifnot(is(metadata, 'data.frame'))

  return(SCData(umi_counts = umi_counts, metadata = metadata))
}

######## -----------------------------------------------------------------------
########                 Private Functions
######## ---------------------------------------------------------------

###### -----------------------------------------
######     .scd_set_count_derived_metadata

.scd_set_count_derived_metadata = function(umi_counts, metadata) {

  if (!is.null(metadata)) {
    # retain metadata columns that are not about to be recreated
    colnames_keep = !(colnames(metadata) %in% c('nGene', 'nUMI', 'mt_prop'))
    metadata = metadata[,colnames_keep,drop=F]
  }

  # create count derived metadata
  nUMI = apply(umi_counts, 2, sum)
  nGene = apply(umi_counts != 0, 2, sum)

  mito_genes = grep(pattern="^mt-", x=rownames(umi_counts), ignore.case=T, value=T)
  mt_prop = apply(umi_counts[mito_genes,,drop=F], 2, sum) / nUMI

  # combine metadata
  if (!is.null(metadata)) {

    metadata = data.frame(nUMI=nUMI, nGene=nGene, mt_prop=mt_prop, metadata)

  } else {

    metadata = data.frame(nUMI=nUMI, nGene=nGene, mt_prop=mt_prop)
  }

  # remove unused levels
  metadata = droplevels(metadata)

  return(metadata)
}

######## -----------------------------------------------------------------------
########                 Public Functions
######## ---------------------------------------------------------------

###### ------------------------------
######     scdata_limit_to_cells

scdata_limit_to_cells = function(scdata, cells_to_keep) {

  # check input types
  stopifnot(is(scdata, 'SCData'))
  stopifnot(is(cells_to_keep, 'character'))

  stopifnot(length(cells_to_keep) > 0)
  stopifnot(all(cells_to_keep %in% colnames(scdata@umi_counts)))

  return(SCData(umi_counts=scdata@umi_counts[,cells_to_keep,drop=F],
                metadata=scdata@metadata[cells_to_keep,,drop=F]))
}

###### ---------------------
######     scdata_limit

scdata_limit = function(scdata, cells_to_keep, genes_to_keep) {

  # check input types
  stopifnot(is(scdata, 'SCData'))
  stopifnot(is(cells_to_keep, 'character'))
  stopifnot(is(genes_to_keep, 'character'))

  stopifnot(length(cells_to_keep) > 0)
  stopifnot(length(genes_to_keep) > 0)
  stopifnot(all(cells_to_keep %in% colnames(scdata@umi_counts)))
  stopifnot(all(genes_to_keep %in% rownames(scdata@umi_counts)))

  return(SCData(umi_counts=scdata@umi_counts[genes_to_keep,cells_to_keep,drop=F],
                metadata=data.frame(scdata@metadata[cells_to_keep,,drop=F])))
}

###### --------------------------
######     scdata_downsample

scdata_downsample = function(scdata, max_umi=NULL) {

  # check input types
  stopifnot(is(scdata, 'SCData'))
  if(!is.null(max_umi)) { stopifnot(is(max_umi, 'numeric')) }

  # determine the max umi to use
  if (is.null(max_umi)) { max_umi = median(scdata@metadata$nUMI) }

  stopifnot(max_umi > 0)
  stopifnot(all(rownames(scdata@metadata)==colnames(scdata@umi_counts)))

  # downsample the umi counts
  umi_counts = Seurat::SampleUMI(scdata@umi_counts, max.umi = floor(max_umi))
  colnames(umi_counts) = colnames(scdata@umi_counts)
  rownames(umi_counts) = rownames(scdata@umi_counts)

  # count derived metadata recalculated upon SCData creation

  return(SCData(umi_counts=umi_counts, metadata=scdata@metadata))
}

###### --------------------------
######     scdata_split

scdata_split = function(scdata, factor) {

  # check input types
  stopifnot(is(scdata, 'SCData'))
  stopifnot(is(factor, 'character'))

  stopifnot(all(rownames(scdata@metadata)==colnames(scdata@umi_counts)))
  stopifnot(factor %in% colnames(scdata@metadata))

  # create an scd for each level of the factor
  scdata_list = list()
  for (v in unique(scdata@metadata[,factor])) {
    v_cell_names = rownames(scdata@metadata)[scdata@metadata[,factor]==v]
    scdata_list[[v]] = scdata_limit_to_cells(scdata = scdata,
                                             cells_to_keep = v_cell_names)
  }

  return(scdata_list)
}

###### ----------------------
######     scdata_combine

scdata_combine = function(scdata_list) {

  # check input types
  stopifnot(is(scdata_list, 'list'))
  for (i in 1:length(scdata_list)) { stopifnot(is(scdata_list[[i]], 'SCData')) }

  # check all gene and metadata names
  for (i in 2:length(scdata_list)) {
    stopifnot(all(rownames(scdata_list[[1]]@umi_counts)==rownames(scdata_list[[i]]@umi_counts)))
    stopifnot(all(colnames(scdata_list[[2]]@metadata)==colnames(scdata_list[[i]]@metadata)))
  }

  # get first scdata's data
  umi_counts = scdata_list[[1]]@umi_counts
  metadata = scdata_list[[1]]@metadata

  # add the other scdatas' data
  scdata_list = scdata_list[-1]
  for (scdata in scdata_list) {
    umi_counts = cbind(umi_counts, scdata@umi_counts)
    metadata = rbind(metadata, scdata@metadata)
  }

  return(SCData(umi_counts=umi_counts, metadata=metadata))
}

###### ------------------------------
######     scdata_save

scdata_save = function(scdata, save_dir) {

  # check input types
  stopifnot(is(scdata, 'SCData'))
  stopifnot(is(save_dir, 'character'))

  # avoid trailing slashes
  save_dir = paste(strsplit(save_dir, split = '/')[[1]], collapse = '/')

  # ensure the save directory exists
  stopifnot(dir.exists(save_dir))

  # save the corrected count matrix
  umi_counts_filename = paste(save_dir, 'umi_counts.csv.bz2', sep = '/')
  if (file.exists(umi_counts_filename)) { stop('umi counts file exists already') }

  write.csv(
    as.matrix(scdata@umi_counts),
    file = bzfile(umi_counts_filename),
    row.names = T
  )

  # save the metadata
  metadata_filename = paste(save_dir, 'metadata.csv', sep = '/')
  if (file.exists(metadata_filename)) { stop('metadata file exists already') }

  write.csv(scdata@metadata,
            file = metadata_filename,
            row.names = T)
}
