


# SEVILLAN PLOT: =========================================================

# rename stuff
library( plyr )
dimnames( pVals_iViC )[[1]] <- mapvalues( dimnames( pVals_iViC )[[1]],
                                          c( "Status", "e4_positive", "White.Matter", "Hippocampus", "Enthorhinal", "ER_Thickness",
                                             "Cer.d18.1.16.0.", "Cer.d18.1.18.0.", "Cer.d18.1.20.0.", "Cer.d18.1.22.0.", "Cer.d18.1.24.0.", "Cer.d18.1.24.1.",
                                             "PC36.5", "PC38.6", "PC40.6",
                                             "Cers1", "Cers2", "PCs" ),
                                          c( "status", "APOE", "white matter", "hippocampus", "enthorhinal", "ER thickness",
                                             "Cer16:0", "Cer18:0", "Cer20:0", "Cer22:0", "Cer24:0", "Cer24:1",
                                             "PC36:5", "PC38:6", "PC40:6",
                                             "light cers.", "heavy cers.", "PCs" ) );
dimnames( pVals_iViC )[[2]] <- mapvalues( dimnames( pVals_iViC )[[2]],
                                          c( "Status", "e4_positive", "White.Matter", "Hippocampus", "Enthorhinal", "ER_Thickness",
                                             "Cer.d18.1.16.0.", "Cer.d18.1.18.0.", "Cer.d18.1.20.0.", "Cer.d18.1.22.0.", "Cer.d18.1.24.0.", "Cer.d18.1.24.1.",
                                             "PC36.5", "PC38.6", "PC40.6",
                                             "Cers1", "Cers2", "PCs" ),
                                          c( "status", "APOE", "white matter", "hippocampus", "enthorhinal", "ER thickness",
                                             "Cer16:0", "Cer18:0", "Cer20:0", "Cer22:0", "Cer24:0", "Cer24:1",
                                             "PC36:5", "PC38:6", "PC40:6",
                                             "light cers.", "heavy cers.", "PCs" ) );

# # Plot results
# source( "FUNCTIONS_utils.R" )
# library( reshape )
# pVals_iViC[ is.na( pVals_iViC ) ] <- 1;
# pVals_iRiC <- melt( pVals_iViC );
# pVals_iRiC$pvalD <- discretePvals( pVals_iRiC[,"value"],
#                                    c( 0.01, 0.001, 0.0001 ) );
# ggplot( pVals_iRiC ) +
#   aes( x = X1,
#        y = X2,
#        color = pvalD ) +
#   geom_point( size = 5 ) +
#   theme( axis.text.x  = element_text( angle=90, 
#                                       hjust=1,
#                                       vjust = 0.5 ) ) +
#   scale_colour_manual( values = c( "white", "yellow", "red", "purple" ) );
# 
# 
# 

# Plot results
library( ggplot2 )
ggplot( p_iXXiC ) +
  aes( x = nXs,
       y = nYs,
       fill = pDics ) +
  geom_raster( ) +
  scale_fill_manual( name = "p value",
                     values = c( hsv( 0, 0.0, 1, 0.5 ),
                                 hsv( 0, 0.4, 1 ),
                                 hsv( 0, 0.6, 1 ),
                                 hsv( 0.8, 1.0, 1 ) ) ) +
  theme( axis.text.x = element_text( angle = 90,
                                     vjust = 0.5,
                                     hjust = 1 ) ) +
  facet_grid( measure ~ . );



ggsave( filename = "metabolicsGLM_v8.pdf",
        width = 5,
        height = 4.2 )


# TRANSLATE STUPID GENES ========================================================


library (org.Hs.eg.db)
library(KEGGgraph)
library(KEGG.db)
library(reactome.db)
pathNames_jGiP <- mget( as.character( allEntrez_iB ), 
                        org.Hs.egPATH,
                        ifnotfound = NA );


