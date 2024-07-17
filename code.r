#1. Species accumulation curve
library(vegan)
library(ggvegan)
# Read the data for the first region
otu1 <- read.delim('ALE_otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- otu1[, -dim(otu1)[2]]
otu1 <- t(otu1)
# Read the data for the second region
otu2 <- read.delim('GKM_otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu2 <- otu2[, -dim(otu2)[2]]
otu2 <- t(otu2)
# Calculate the species accumulation curves
sp1 <- specaccum(otu1, method = 'random')
sp2 <- specaccum(otu2, method = 'random')
# Create a new plot
plot(sp1, ci.type = 'poly', col = '#009e73', lwd = 2, ci.lty = 0, ci.col = FALSE)
boxplot(sp1, col = 'yellow', add = TRUE, pch = '+')
#poolaccum() uses some popular methods to estimate the number of these unseen species and adds them to the observed species richness
#details ?poolaccum ??autoplot
pool1 <- poolaccum(otu1)
summary(pool1, display = 'chao')
autoplot(pool1)
autoplot(pool1, facet = FALSE)
# Add the curves from the second region to the existing plot
lines(sp2, ci.type = 'poly', col = '#e69f00', lwd = 2, ci.lty = 0, ci.col = FALSE)
boxplot(sp2, col = '#e69f00', add = TRUE, pch = '+')

pool2 <- poolaccum(otu2)
summary(pool2, display = 'chao')
autoplot(pool2)
autoplot(pool2, facet = FALSE)

#2. Calculation of alpha diversity
library(vegan)
library(picante)

# Read in the OTU table
otu <- read.delim('ALE_otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
# Transpose the OTU data
otu <- t(otu)

alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon', base = 2)
  Simpson <- diversity(x, index = 'simpson')    # Note: This is the Gini-Simpson index
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  # Keep four decimal places
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  
  result <- data.frame(observed_species, ACE, Chao1, Shannon, Simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE, Chao1, Shannon, Simpson,
                         PD_whole_tree, goods_Coverage)
  }
  
  result
}

# No need to calculate phylogenetic diversity
alpha <- alpha_diversity(otu)
# Output the result and save locally
write.csv(alpha, 'alpha_diversity.csv', quote = FALSE)


#3. Regression analysis and mantel test
library(vegan)

# Read the matrix directly
dist.abund <- as.dist(read.delim('jaccard.txt', row.names = 1, sep = '\t', check.names = FALSE))
df <- read.csv('environment.csv', header= TRUE)

# Use environmental variables and standardize the data to eliminate differences in scales
env <- df[ ,2]
scale.env <- scale(env, center = TRUE, scale = TRUE)
dist.env <- dist(scale.env, method = 'euclidean')

# Calculate the actual geographic distance between plots based on longitude and latitude
library(geosphere)
geo <- data.frame(df$Longitude, df$Latitude)
d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)

# Correlate species matrix with environmental variables using Spearman correlation and 9999 permutations for significance testing
abund_env <- mantel(dist.abund, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_env

# Combine distances into a dataframe
aa <- as.vector(dist.abund)
ss <- as.vector(dist.env)
gg <- as.vector(dist.geo)
mydata <- data.frame(aa, ss)

# Fit a linear model and summarize the fit
myfit <- lm(dist.abund ~ dist.env, data = mydata)
summary(myfit)

# Extract important values
coefficients(myfit)[1]  # Intercept
coefficients(myfit)[2]  # Slope of site_dis_km
summary(myfit)$adj.r.squared  # Adjusted R2


# Perform multiple regression on distance matrices
library(ecodist)

mrm_fit <- MRM(dist.abund ~ dist.env, nperm = 999)

# Summarize the fit
summary(mrm_fit)

# Extract important values
mrm_fit$coef  # Coefficients
mrm_fit$r.squared  # R-squared


# Scatter plot of the correlation between species abundance distance and environmental distance with color representing geographic distance
library(ggplot2)
mm <- ggplot(mydata, aes(y = aa, x = ss)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
  geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, aes(fill = gg / 1000)) + 
  labs(x = "Environmental Distance", y = "jaccard", fill = "Physical Separation (km)") + 
  theme(axis.text.x = element_text(face = "bold", colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "top",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "#e69f00", low = "#009e73")
mm

#4. UPGMA

library(ggtree)
library(treeio)
library(ape)

# Load environmental data
df <- read.csv('environment.csv', header = TRUE)

# Use environmental data and standardize them to eliminate dimensional differences
env <- df[, 2:4]
scale.env <- scale(env, center = TRUE, scale = TRUE)
dist.env <- dist(scale.env, method = 'euclidean')

# Load taxonomic abundance table
taxo <- read.delim('18s_order.txt', row.names = 1, sep = '\t')

# Perform hierarchical clustering based on the distance matrix
tree <- hclust(dist.env, method = 'average')

# Convert the clustering result to phylogenetic format
df_tree <- as.phylo(tree)
write.tree(df_tree, file = "tree.nwk")

#5. Sloan's neutral model

library(Hmisc)
library(minpack.lm)
library(stats4)

# Read species or taxon abundance table, where rows are taxa and columns are samples
spp <- read.csv('abundance.txt', head = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
spp <- t(spp)

# Fit Sloan's neutral model (2006) to the abundance table and return fit statistics
# Fit model parameters using Non-linear Least Squares (NLS)
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m / N
spp.bi <- 1 * (spp > 0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by = 0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1 / N
m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
m.fit  # Get m value
m.ci <- confint(m.fit, 'm', level = 0.95)
freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
Rsqr  # Get R2 value of the model

# Output three statistical result tables: mean relative abundance (p.csv), observed occurrence frequency (freq.csv), and predicted occurrence frequency (freq.pred.csv)
write.csv(p, file = "p.csv")
write.csv(freq, file = "freq.csv")
write.csv(freq.pred, file = "freq.pred.csv")

# p is the mean relative abundance
# freq is the observed occurrence frequency
# freq.pred is the predicted occurrence frequency, i.e., the fitted value from the neutral model

# Plot the results
bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[,2:3])
inter.col <- rep('black', nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#A52A2A'  # Frequency below neutral community model prediction
inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#29A6A6'  # Frequency above neutral community model prediction

library(grid)
grid.newpage()
pushViewport(viewport(h = 0.6, w = 0.6))
pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), yData = c(0, 1.02), extension = c(0.02, 0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = 'blue', lwd = 2), default = 'native')
grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = 'blue', lwd = 2, lty = 2), default = 'native')
grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = 'blue', lwd = 2, lty = 2), default = 'native')
grid.text(y = unit(0, 'npc') - unit(2.5, 'lines'), label = 'Mean Relative Abundance (log10)', gp = gpar(fontface = 2))
grid.text(x = unit(0, 'npc') - unit(3, 'lines'), label = 'Frequency of Occurrence', gp = gpar(fontface = 2), rot = 90)

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr =", round(Rsqr, 3), "\n", "Nm =", round(coef(m.fit) * N)), x = x[j], y = y[i], just = just)
}
x <- unit(1:4 / 5, "npc")
y <- unit(1:4 / 5, "npc")
draw.text(c("centre", "bottom"), 4, 1)
