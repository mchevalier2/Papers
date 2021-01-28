library(crestr)
library(pals)
library(viridis)
library(raster)
library(sp)

{
    EXT=c(10,42,-40,-12)
    data(M1)
    M1 <- sp::spTransform(M1, CRS("+proj=longlat +datum=WGS84"))
    M1=raster::crop(M1,EXT)

    variables=c("Tmean_Wet_Q","Tmin_Cold_M","Prec_Warm_Q","Prec_Cold_Q","Aridity")


    VARIABLE.COL=list("Tmin_Cold_M"=pals::coolwarm(27),
                      "Tmean_Wet_Q"=pals::coolwarm(21),
                      "Aridity"=pals::ocean.speed(22),
                      "Prec_Warm_Q"=rev(pals::ocean.delta(30)[1:15]),
                      "Prec_Cold_Q"=rev(pals::ocean.delta(34)[1:17]),
                      "distrib"=rev(viridis::viridis(7))
                    )

    VARIABLE.COL1000=list("Tmin_Cold_M"=pals::coolwarm(500),
                          "Tmean_Wet_Q"=pals::coolwarm(500),
                          "Aridity"=pals::ocean.speed(500),
                          "Prec_Warm_Q"=rev(pals::ocean.delta(1000)[1:500]),
                          "Prec_Cold_Q"=rev(pals::ocean.delta(1000)[1:500])
                        )

    COL.BIOMES=list("Deserts and xeric shrublands"=rgb(230,171,2,maxColorValue=255),
                    'Mediterranean Forests, woodlands and scrubs'=rgb(117,112,179,maxColorValue=255),
                    'Montane grasslands and shrublands'=rgb(231,41,138,maxColorValue=255),
                    'Flooded grasslands and savannas'=rgb(166,118,29,maxColorValue=255),
                    'Tropical and subtropical grasslands, savannas and shrublands'=rgb(217,95,2,maxColorValue=255),
                    'Tropical and subtropical moist broadleaf forests'=rgb(27,158,119,maxColorValue=255),
                    'Tropical and subtropical dry broadleaf forests'=rgb(102,166,30,maxColorValue=255),
                    'Mangroves'=rgb(252,205,229, maxColorValue=255),
                    'Lakes'='black'
                  )

    CLASS_WIDTH=list("Prec_Warm_Q"=50,
                     "Prec_Cold_Q"=20,
                     "Tmean_Wet_Q"=1,
                     "Aridity"=0.05,
                     'Tmin_Cold_M'=1
                    )

    WC_NAMES=list("Prec_Warm_Q"='bio18',
                  "Prec_Cold_Q"='bio19',
                  "Tmean_Wet_Q"='bio8',
                  "Aridity"='ai',
                  'Tmin_Cold_M'='bio6'
                 )

    VARIABLES=list()
    for(v in variables)  {
        VARIABLES[[v]]=dbRequest(paste0( "SELECT DISTINCT longitude, latitude, ",WC_NAMES[[v]],
                                         " FROM wc_qdgc",
                                         " WHERE latitude >= ", EXT[3],
                                           " AND latitude <= ", EXT[4],
                                           " AND longitude >= ", EXT[1],
                                           " AND longitude <= ", EXT[2],
                                           " AND ai IS NOT NULL",
                                           ";"
                                        )
                                )
    }

    VARIABLES[["Aridity"]][,3]=sqrt(VARIABLES[["Aridity"]][,3]/10000)

    VARIABLES_NAMES=list("Prec_Warm_Q"=c("Precipitation of the Warmest Quarter","R","lognormal","PWarmQ","(mm)", "Precipitation\nof the Warmest Quarter (mm)"),
                         "Prec_Cold_Q"=c("Precipitation of the Coldest Quarter","R","lognormal","PColdQ","(mm)", "Precipitation\nof the Coldest Quarter (mm)"),
                         "Aridity"=c("Aridity Index","R","lognormal","Aridity","", "Aridity Index"),
                         "Tmean_Wet_Q"=c("Temperature of the Wettest Quarter","L","normal","TWetQ","(°C)","Temperature\nof the Wettest Quarter (°C)"),
                         "Tmin_Cold_M"=c("Minimum temperature of the coldest month","R","normal","TColdM","(°C)","Minimum temperature\nof the coldest month (°C)")
                       )

    XRANGE=list("Prec_Warm_Q"=c(0,750),
                "Prec_Cold_Q"=c(0,349),
                "Aridity"=c(0,1.14),
                "Tmean_Wet_Q"=c(7,29),
                "Tmin_Cold_M"=c(-6,20)
              )

    RASTERS=list()
    CLIMATE_SPACE=VARIABLES[["Aridity"]][,c(1,2)]
    for(v in variables){
        CLIMATE_SPACE=cbind(CLIMATE_SPACE,VARIABLES[[v]][,3])
        RASTERS[[v]] = VARIABLES[[v]]
        RASTERS[[v]][,3] = ((VARIABLES[[v]][,3]-min(VARIABLES[[v]][,3],na.rm=TRUE))/diff(XRANGE[[v]])*diff(XRANGE[[v]])/CLASS_WIDTH[[v]])%/%1+1
        RASTERS[[v]] = rasterFromXYZ(RASTERS[[v]])
    }

    POLTYPES=POLTYPES_UNIQUE=BIOMES=BIOMES_UNIQUE=TAXONOMY=STATS=list()
    POLLEN_TAXA = unique(rio::import('pollenTypes_AtlasSA.xlsx', skip=1)[,'ProxyName'])
    POLLEN_TAXA = POLLEN_TAXA[!POLLEN_TAXA %in% c('Lactucoideae', 'Mimosoideae')]
    POLLEN_TAXA = POLLEN_TAXA[!POLLEN_TAXA %in% c('Bruguiera', 'Blaeria-type', 'Rhaphiostylis', 'Carpacoce')]
    print('I do not have a good classification for Lactucoideae, Mimosoideae therefore I exclude them for now.')
    print('I do not have a enough points for Bruguiera, Rhaphiostylis and Blaeria therefore I exclude them for now.')

    POLLEN_TAXA = sort(POLLEN_TAXA)

    taxonID2proxy <- data.frame("taxonID" = NA, "proxyName" = NA)
    PSE = rio::import('pollenTypes_AtlasSA.xlsx')
    for (taxLevel in 1:3) {
        for (tax in PSE$ProxyName[ PSE$Level == taxLevel ]) {
            for (w in which(PSE$ProxyName == tax)) {
                taxonIDs <- getTaxonID(
                  PSE$Family[w],
                  PSE$Genus[w],
                  PSE$Species[w],
                  taxaType = 1
                )
                if (length(taxonIDs) > 0) {
                    existingTaxa <- taxonIDs %in% taxonID2proxy[, "taxonID"]
                    # If the taxon was first assigned to higher group, Reassign.
                    if (sum(existingTaxa) > 0) {
                        taxonID2proxy[taxonID2proxy[, "taxonID"] %in% taxonIDs, "proxyName"] <- tax
                    }
                    if (sum(existingTaxa) != length(taxonIDs)) {
                        taxonID2proxy <- rbind(
                          taxonID2proxy,
                          data.frame("taxonID" = taxonIDs[!existingTaxa],
                                     "proxyName" = rep(tax, sum(!existingTaxa)),
                                     stringsAsFactors = FALSE)
                                   )
                    }
                }
            }
        }
    }
    taxonID2proxy = taxonID2proxy[-1, ]


    LIST_OF_TAXONIDS = list()
    for(pol in POLLEN_TAXA){
        LIST_OF_TAXONIDS[[pol]] = paste0('(', paste(as.character(taxonID2proxy[taxonID2proxy[,2] == pol, 1]), collapse=', '),')')
        vv=paste(sapply(variables, function(x) return(WC_NAMES[[x]])), collapse=',')
        dat = dbRequest(paste0( "SELECT distrib_qdgc.longitude,distrib_qdgc.latitude,species,genus,",vv,
                                "  FROM taxa,distrib_qdgc,wc_qdgc ",
                                " WHERE taxa.taxonid=distrib_qdgc.taxonid",
                                "   AND distrib_qdgc.longitude=wc_qdgc.longitude",
                                "   AND distrib_qdgc.latitude=wc_qdgc.latitude",
                                "   AND taxa.taxonID IN ", LIST_OF_TAXONIDS[[pol]],
                                "   AND distrib_qdgc.latitude >= ", EXT[3],
                                "   AND distrib_qdgc.latitude <= ", EXT[4],
                                "   AND distrib_qdgc.longitude >= ", EXT[1],
                                "   AND distrib_qdgc.longitude <= ", EXT[2],
                                "   AND ai IS NOT NULL"
                              )
                            )

        if(nrow(dat) > 0){
            POLTYPES[[pol]]=dat
            colnames(POLTYPES[[pol]])[5:9] = variables
            POLTYPES[[pol]][,"Aridity"]= sqrt(POLTYPES[[pol]][,"Aridity"]/10000)

            BIOMES[[pol]] = dbRequest(paste0("SELECT distrib_qdgc.longitude,distrib_qdgc.latitude,taxa.taxonid,biome",
                                              " FROM taxa,distrib_qdgc,wwf_qdgc ",
                                              "WHERE taxa.taxonid=distrib_qdgc.taxonid",
                                              "  AND distrib_qdgc.longitude=wwf_qdgc.longitude",
                                              "  AND distrib_qdgc.latitude=wwf_qdgc.latitude",
                                              "  AND taxa.taxonID IN ", LIST_OF_TAXONIDS[[pol]],
                                              "  AND distrib_qdgc.latitude >= ", EXT[3],
                                              "  AND distrib_qdgc.latitude <= ", EXT[4],
                                              "  AND distrib_qdgc.longitude >= ", EXT[1],
                                              "  AND distrib_qdgc.longitude <= ", EXT[2]
                                            )
                                          )
            BIOMES_UNIQUE[[pol]]=table(BIOMES[[pol]][,4])
            BIOMES_UNIQUE[[pol]] = round(100 * BIOMES_UNIQUE[[pol]] / sum(BIOMES_UNIQUE[[pol]]),1)
            BIOMES[[pol]]=unique(BIOMES[[pol]][,c(1,2,4)])
            POLTYPES[[pol]] = merge(POLTYPES[[pol]], BIOMES[[pol]], by=c('longitude', 'latitude'))
            POLTYPES_UNIQUE[[pol]]=unique(POLTYPES[[pol]][,c(1,2,5:10)])

            sp = paste(unique(POLTYPES[[pol]]$species), collapse="', '")
            TAXONOMY[[pol]] = dbRequest(paste0("   SELECT kingdom, phylum,class_name, order_name, family,genus,species ",
                                               "     FROM taxa ",
                                               "    WHERE taxa.taxonID IN ", LIST_OF_TAXONIDS[[pol]],
                                               "      AND species IN ('", sp ,"')",
                                               " ORDER BY genus,species"
                                              )
                                            )

            STATS[[pol]]=list()
            for(v in variables)  STATS[[pol]][[v]]=fitPDFS(v, pol)
        } else {
            print(pol)
        }
    }

    VARIABLES[['distrib']] = dbRequest(paste0("SELECT distrib_qdgc.longitude,distrib_qdgc.latitude,count(*), countryname ",
                                              "  FROM taxa,distrib_qdgc,wc_qdgc, geo_qdgc ",
                                              " WHERE taxa.taxonid=distrib_qdgc.taxonid",
                                              "   AND distrib_qdgc.longitude=wc_qdgc.longitude",
                                              "   AND distrib_qdgc.latitude=wc_qdgc.latitude",
                                              "   AND geo_qdgc.longitude=wc_qdgc.longitude",
                                              "   AND geo_qdgc.latitude=wc_qdgc.latitude",
                                              "   AND taxa.taxonID IN ", LIST_OF_TAXONIDS[[pol]],
                                              "   AND distrib_qdgc.latitude >= ", EXT[3],
                                              "   AND distrib_qdgc.latitude <= ", EXT[4],
                                              "   AND distrib_qdgc.longitude >= ", EXT[1],
                                              "   AND distrib_qdgc.longitude <= ", EXT[2],
                                              "   AND ai IS NOT NULL",
                                              " GROUP BY distrib_qdgc.longitude,distrib_qdgc.latitude,countryname"
                                              )
                                            )
    RASTERS[['distrib']] = VARIABLES[['distrib']]
    RASTERS[['distrib']][,3] = log10(RASTERS[['distrib']][,3])
    RASTERS[['distrib']] = rasterFromXYZ(RASTERS[['distrib']][, -4])

    BIOMES2NBs=as.data.frame(matrix(1:9, ncol=1))
    rownames(BIOMES2NBs) = names(COL.BIOMES)
    VARIABLES[['biomes']] = dbRequest(paste0("SELECT DISTINCT wwf_qdgc.longitude,wwf_qdgc.latitude,biome ",
                                             "  FROM wc_qdgc, wwf_qdgc",
                                             " WHERE wwf_qdgc.longitude=wc_qdgc.longitude",
                                             "   AND wwf_qdgc.latitude=wc_qdgc.latitude",
                                             "   AND wwf_qdgc.latitude >= ", EXT[3],
                                             "   AND wwf_qdgc.latitude <= ", EXT[4],
                                             "   AND wwf_qdgc.longitude >= ", EXT[1],
                                             "   AND wwf_qdgc.longitude <= ", EXT[2],
                                             "   AND ai IS NOT NULL"
                                            )
                                  )
    RASTERS[['biomes']] = VARIABLES[['biomes']]
    RASTERS[['biomes']][,3] = BIOMES2NBs[RASTERS[['biomes']][,3],]
    RASTERS[['biomes']] = rasterFromXYZ(RASTERS[['biomes']])



    ROUNDS=list("Tmin_Cold_M"=1,
                "Tmean_Wet_Q"=1,
                "Aridity"=3,
                "Prec_Warm_Q"=0,
                "Prec_Cold_Q"=0
              )
}
