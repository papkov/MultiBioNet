alz_intact <- read.csv('../data/alz1_features.csv')

alz.hex <- as.h2o(alz_intact, destination_frame = 'train.hex')

alz.dl <- h2o.deeplearning(x = colnames(alz_intact)[-c(1, 2)], 
                           training_frame = alz.hex,
                           autoencoder = T,
                           reproducible = T,
                           seed = 992,
                           hidden = c(6,5,6), epochs = 50)

alz.anom <- h2o.anomaly(alz.dl, alz.hex, per_feature = F)
head(alz.anom)
err <- as.data.frame(alz.anom)

plot(sort(err$Reconstruction.MSE), main='Reconstruction Error')

alz_intact_auto <- alz_intact[err$Reconstruction.MSE < 0.1,]

write.csv(alz_intact_auto, '../data/alz_autoenc_curated_features.csv')
