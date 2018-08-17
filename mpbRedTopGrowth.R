
# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "mpbRedTopGrowth",
  description = "Mountain Pine Beetle Red Top Growth Model: Short-run Potential for Establishment, Eruption, and Spread",
  keywords = c("mountain pine beetle, outbreak dynamics, eruptive potential, spread, climate change, twitch response"),
  authors = c(
    person(c("Barry", "J"), "Cooke", email = "barry.cooke@canada.ca", role = c("aut", "cre")),
    person(c("Alex", "M"), "Chubaty", email = "alexander.chubaty@canada.ca", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  reqdPkgs = list("amc", "data.table", "ggplot2", "quickPlot", "raster", "reproducible"),
  parameters = rbind(
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the interval between save events"),
    defineParameter(".useCache", "numeric", FALSE, NA, NA, "Should this entire module be run with caching activated?"),
    defineParameter("dataset", "character", "Boone2001", NA, NA, "Which dataset to use for stand dynamic model fitting. One of 'Boone2001' (default), 'Berryman1979_fit', or 'Berryman1979_forced'. Others to be implemented later."),
    defineParameter("growthInterval", "numeric", 1, NA, NA, "This describes the interval time between growth events")
  ),
  inputObjects = bind_rows(
    expectsInput("climateSuitabilityMap", "RasterLayer", "A climatic suitablity map for the current year."),
    expectsInput("massAttacksDT", "data.table", "Current MPB attack map (number of red attacked trees)."),
    expectsInput("massAttacksMap", "RasterStack", "Historical MPB attack maps (number of red attacked trees)."),
    expectsInput("pineDT", "data.table", "Current lodgepole and jack pine available for MPB."),
    expectsInput("pineMap", "data.table", "Current lodgepole and jack pine available for MPB.")
  ),
  outputObjects = bind_rows(
    createsOutput("massAttacksDT", "data.table", "Current MPB attack map (number of red attacked trees)."),
    createsOutput("pineDT", "data.table", "Current lodgepole and jack pine available for MPB.")
  )
))

## event types
#   - type `init` is required for initiliazation

doEvent.mpbRedTopGrowth <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
    "init" = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- sim$mpbRedTopGrowthInit(sim)
      sim <- sim$mpbRedTopGrowthPlotInit(sim)

      # schedule future event(s)
      #sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopGrowth", "plot")
    },
    "grow" = {
      # do stuff for this event
      sim <- sim$mpbRedTopGrowthGrow(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopGrowth", "grow")
    },
    "plot" = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      sim <- sim$mpbRedTopGrowthPlot(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopGrowth", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    "save" = {
      rtmp <- update(rtmp, cell = sim$massAttacks[, ID], v = sim$massAttacks[, RedTrees])
      writeRaster(r, filename = file.path(outputPath(sim), paste0("massAttacks", time(sim), ".tif")))
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  if (!('studyArea' %in% sim$.userSuppliedObjNames)) {
    f <- file.path(modulePath(sim), "mpbRedTopGrowth", "data", "studyArea.kml")
    prj <- paste("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113",
                 "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
    sim$studyArea <- readOGR(f, "studyArea.kml") %>%
      sp::spTransform(., prj)
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

mpbRedTopGrowthInit <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  ## create a data.table consisting of the reduced map of current MPB distribution,
  ## proportion pine, and climatic suitability;
  ## use only the start year's non-zero and non-NA data
  r <- sim$massAttacksMap[[paste0("X", start(sim))]]
  ids <- which(!is.na(r[]) | (r[] > 0))
  mpb.sp <- raster::xyFromCell(r, cell = ids)
  sim$massAttacksDT <- data.table(
    ID = ids,
    #X = mpb.sp[, 1],
    #Y = mpb.sp[, 2],
    NUMTREES = r[ids],
    CLIMATE = raster::extract(sim$mpbClimateDataMaps, mpb.sp)
  )
  sim$massAttacksDT[NUMTREES > 0]
  rm(r)

  ## growth data
  sim$growthData <- switch(P(sim)$dataset,
    "Berryman1979_fit" = {
      ## Berryman1979_forced
      data.frame(
        year = c(1:15),
        log10Xtm1 = c(-3.1, -2.75, -2.7, -2.4, -2.3, -1.2, -1, 0.2, 0.9, 0.65,
                      1.05, 0.95, 1.1, 1.5, 1.85),
        log10Rt = c(0.35, 0.4, 0.1, -0.4, -0.65, 0.3, 1, 0.75, 1.2, -0.7,
                    -0.4, 0.2, 0.45, 0.3, -0.78),
        study = c(rep("Tunnock 1970", 9), rep("Parker 1973", 6)),
        stringsAsFactors = TRUE
      )
    },
    "Berryman1979_forced" = {
      ## same as Berryman1979_fit
      data.frame(
        year = c(1:15),
        log10Xtm1 = c(-3.1, -2.75, -2.7, -2.4, -2.3, -1.2, -1, 0.2, 0.9, 0.65,
                    1.05, 0.95, 1.1, 1.5, 1.85),
        log10Rt = c(0.35, 0.4, 0.1, -0.4, -0.65, 0.3, 1, 0.75, 1.2, -0.7,
                      -0.4, 0.2, 0.45, 0.3, -0.78),
        study = c(rep("Tunnock 1970", 9), rep("Parker 1973", 6)),
        stringsAsFactors = TRUE
      )
    },
    "Boone2001" = {
      data <- read.csv(file.path(modulePath(sim), "mpbRedTopGrowth", "data", "BooneCurveData2.csv"))
      data$Site <- c(rep("A", 6), rep("B", 6), rep("D", 5), rep("E", 4), rep("F", 4), rep("G", 3))
      data$Year <- c(2000:2005, 2000:2005, 2001:2005, 2002:2005, 2002:2005, 2003:2005)
      data
    }
  )

  ## define growth function (from regression) for each dataset
  sim$growthFunction <- switch(P(sim)$dataset,
     "Berryman1979_fit" = {
       function(x, s) {
         # TODO: check this works
         m <- lm(log10Rt ~ poly(log10Xtm1, 3, raw = TRUE), data = sim$growthData)
         s * unname(predict(m, newdata = data.frame(log10Xtm1 = x)))
       }
     },
     "Berryman1979_forced" = {
       function(x, s) {
         # TODO: check this works
         poly3.params <- c(1.1, -0.2, -0.9, -0.24)
         s * (poly3.params[4] * x^3 + poly3.params[3] * x^2 + poly3.params[2] * x + poly3.params[1])
       }
     },
     "Boone2001" = {
       function(x, s) {
         ## x is number of attacked trees (NUMTREES)
         ## s is scaling parameter (0,1), based on climate (CLIMATE)

         ## mortality from emigration/dispersal
         # r: relative stocking value (0,1)
         # d: slope parameter [1,Inf)
         # s: scaling parameter (0,1)
         m_e <- function(r, d, s) {
           s * exp(1 - d * r)
         }

         # use 2004 data as baseline for unweakened hosts (i.e., a good year for trees)
         m <- lm(amc::logit(PropKilled) ~ log(Attacked), data = subset(sim$growthData, Year == "2004"))
         a <- 0.9              # scale parameter; TODO: explain this
         d <- 3                # d: slope parameter [1,Inf)
         r <- 0.2              # r: relative stocking value (0,1)
         fudge2 <- 0.9         # from MacQuarrie 2011 (Fig 3d); TODO: extract from raw data
         fudge <- fudge2 + 0.3 # somewhat arbitrary; chosen so that the resulting curve passes 1 when flexed

         # resulting function
         log(amc::hill(m$coefficients[[1]], m$coefficients[[2]], exp(a * x))) +
           (fudge - m_e(r, d, s) - 0.03 * exp(a * x))
         }
     }
  )
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

mpbRedTopGrowthPlotInit <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event

  ### see ggplot docs at http://docs.ggplot2.org/current/
  gg <- if (grepl("Berryman1979", P(sim)$dataset)) {
    ggplot(sim$growthData) +
      geom_point(aes(x = log10Xtm1, y = log10Rt, shape = study)) +
      scale_shape(solid = FALSE) +
      xlim(-3.2, 2) + ylim(-1.5, 1.5) +
      labs(title = switch(P(sim)$dataset,
                          "Berryman1979_fit" = "Berryman (1979) [fit]",
                          "Berryman1979_forced" = "Berryman (1979) [forced]"
                  ),
           x = "X[t-1] (log10 trees/ha/yr)",
           y = "R[t] = log10 x[t]/x[t-1]") +
      geom_hline(aes(yintercept = 0)) +
      switch(P(sim)$dataset,
        "Berryman1979_fit" = {
          stat_smooth(aes(x = log10Xtm1, y = log10Rt), method = "lm",
                      formula = y ~ poly(x, 3, raw = TRUE))
        },
        "Berryman1979_forced" = {
          stat_function(fun = sim$growthFunction, colour = "blue")
        }
      )
  } else if (P(sim)$dataset == "Boone2001") {
    ggplot(sim$growthData) +
      xlim(-3.2, 6) + ylim(-1.5, 1.5) +
      labs(title = "Boone et al. (2001)",
           x = "X[t-1] (log trees/ha/yr)",
           y = "R[t] = log x[t]/x[t-1]") +
      geom_hline(aes(yintercept = 0)) +
      stat_function(fun = sim$growthFunction, args = list(s = 0.9), colour = "purple")
  }

  ### save the object to the simList
  sim$mpbRedTopGrowthPlotGG <- gg

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

mpbRedTopGrowthPlot <- function(sim) {
  currentAttack <- amc::dt2raster(sim$massAttacksDT, sim$massAttacksMap, "NUMTREES")
  Plot(currentAttack, addTo = "sim$massAttacksMap")

  currentPine <- amc::dt2raster(sim$massAttacksDT, sim$massAttacksMap, "PROPPINE")
  Plot(currentPine, addTo = "sim$massAttacksMap")

  scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "mpbRedTopGrowth", "plot")

  return(invisible(sim))
}

mpbRedTopGrowthGrow <- function(sim) {
  ## determine the actual growth based on the actual number of attacked trees/ha
  xt <- function(xtminus1, cs) {
    map.res <- xres(xtminus1)
    per.ha <- 10^sim$growthFunction(log10(xtminus1), cs) * xtminus1 ## TODO: something is off about this
    return(map.res * per.ha)
  }

  sim$massAttacksDT <- sim$massAttacksDT[NUMTREES := xt(NUMTREES, CLIMATE)]
}
