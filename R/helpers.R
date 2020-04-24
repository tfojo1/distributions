

##-----------------------------##
##-- Matching Variable Names --##
##-----------------------------##


match.variables <- function(dist, values, values.error.name = '')
{
    was.vector = is.null(dim(values))
    if (was.vector)
    {
        if (dist@n.var==1 && is.null(names(values)))
            values = as.matrix(values)
        else
            values = t(as.matrix(values))
    }

    if (is.null(dist@var.names) || is.null(dimnames(values)[[2]]))
    {
        if (dist@n.var==1)
        {}
        else if (dim(values)[2] != dist@n.var)
        {
            if (is.null(dist@var.names))
            {
                if (was.vector)
                    stop(paste0("Since variables are not named for the distribution, ",
                         "there must be ",
                         ifelse(dist@n.var==1, "one value", paste0(dist@n.var, " values")),
                         ifelse(values.error.name=='', " passed", paste0("in '", values.error.name, "'"))
                         ))
                else
                    stop(paste0("Since variables are not named for the distribution, ",
                                "there must be ",
                                ifelse(dist@n.var==1, "one cplumn", paste0(dist@n.var, " columns")),
                                ifelse(values.error.name=='', "in the given values", paste0("in '", values.error.name, "'"))
                    ))
            }
            else
            {
                if (vas.vector)
                    stop(paste0(ifelse(values.error.name=='',
                                       "Since the given values are not named, they must be of length ",
                                       paste0("Since '", values.error.name, "' does not have named values, '", values.error.name, "' must be of length ")),
                                dist@n.var
                                ))
                else
                    stop(paste0(ifelse(values.error.name=='',
                                       "Since the columns in the given values are not named, there must be exactly ",
                                       paste0("Since '", values.error.name, "' does not have named columns, '", values.error.name, "' must have exactly ")),
                                dist@n.var, " columns"
                    ))
            }

        }

        dimnames(values) = list(NULL, variable=dist@var.names)

        values
    }
    else
    {
        missing.variables = setdiff(dist@var.names, dimnames(values)[[2]])
        if (length(missing.variables)>1)
            stop(paste0("Values for the following variable",
                        ifelse(length(missing.variables)==1, ' is', 's are'),
                        " not present",
                        ifelse(values.error.name=='', '', paste0(" in '", values.error.name, "' ")),
                        ": ",
                        paste0("'", missing.variables, "'", collapse=', ')))

        matrix(values[,dist@var.names], ncol=dist@n.var, dimnames=list(NULL, variable=dist@var.names))
    }
}

match.values.to.variables <- function(values, var.names, n.var=length(var.names))
{

}

map.keep.indices <- function(dist, vars.to.keep)
{
    if (class(vars.to.keep)=='character')
    {
        if (is.null(dist@var.names))
            stop("Cannot subset a distribution by variable names in vars.to.keep when variables are not named")

        vars.to.keep = unique(vars.to.keep)
        erroneous.variables = setdiff(vars.to.keep, dist@var.names)
        if (length(erroneous.variables)>1)
            stop(paste0("The following ",
                        ifelse(length(erroneous.variables)==1,
                               'is not the name of a variable',
                               'are not names of variables'),
                        " in the distribution: ",
                        paste0("'", erroneous.variables, "'", collapse=', ')))

        indices = 1:dist@n.var
        names(indices) = dist@var.names
        sort(indices[vars.to.keep])
    }
    else if (class(vars.to.keep)=='logical')
    {
        if (length(vars.to.keep) != dist@n.var)
            stop(paste0("If 'vars.to.keep is a logical vector, it must be of length ",
                        dist@n.var, ", with a value for each variable"))

        if (!any(vars.to.keep))
            stop("In subsetting a distribution, at least one element of 'vars.to.keep' must be true")

        (1:dist@n.var)[vars.to.keep]
    }
    else if (class(vars.to.keep)=='integer' || class(vars.to.keep)=='numeric')
    {
        if (any(round(vars.to.keep)!=vars.to.keep))
            stop("If using numeric indices, all elements of 'vars.to.keep' must be integers")

        sort(unique(vars.to.keep))
    }
    else
        stop(paste0("'vars.to.keep' must be either a character vector of variable names, an integer or numeric vector of indices representing variables, or a logical vector of length ",
                    dist@n.var, " indicating which variables to keep"))
}

##----------------------------------------------##
##-- Getting Inetrvals/Quantiles From Samples --##
##----------------------------------------------##

do.get.highest.density.interval.from.samples <- function(samples, coverage, weights=1)
{
    stop('Highest-density intervals have not been implemented')
}

##-----------##
##-- OTHER --##
##-----------##

print.improper.warnings <- function(dist, description="generate random samples from")
{

    if (any(dist@is.improper))
    {
        warning(paste0("Cannot ", description, " improper variable",
                       ifelse(sum(dist@is.improper)==1, ' ', 's '),
                       ifelse(!is.null(dist@var.names),
                              paste0("'", dist@var.names[dist@is.improper], "'", collapse=", "),
                              ''),
                       ": NAs will be generated"))
    }
}

log.sum.exp <- function(components)
{
    max.component = max(components)
    max.component + log(sum(exp(components - max.component)))
}

is.supported.by.dist <- function(dist, x, aggregate.variables=F)
{
    rv = sapply(1:dist@n.var, function(i){
        is.supported(dist@support[[i]], x[,i])
    })

    if (aggregate.variables)
        apply(rv, 1, all)
    else
        rv
}
