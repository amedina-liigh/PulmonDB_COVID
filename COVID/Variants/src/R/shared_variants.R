library(tidyverse)
library(ComplexHeatmap)
library(UpSetR)

# Read data
df <- read_delim(here::here("COVID/Variants/data/variantes-COVID.txt"))

# List of unique SNPs
df.unique <- pivot_longer(df, cols = 1:4) %>% 
  filter(!duplicated(value) & !is.na(value)) %>% 
  rename( Study = name, Variant = value)

write.csv(df.unique, 
          here::here("COVID/Variants/output/unique_variants_COVID.csv"), 
          quote = F,
          row.names = F
          )

# Unique SNPs
lt <- apply(df, 2, unique)
lt <- lapply(lt, na.omit)

m <- make_comb_mat(lt)

# Upset

pdf(here::here("COVID/Variants/fig/UpSet_variants_COVID.pdf"),
    width = 8,
    height = 3
    )
ht <- draw(UpSet(m, pt_size = unit(5, "mm"), lwd = 3,
                 comb_col = c("red", "black", "blue")[comb_degree(m)],
                 right_annotation = rowAnnotation(
                   "Set size" = anno_barplot(set_size(m), 
                                             border = FALSE, 
                                             #height = unit(2, "cm"),
                                             gp = gpar(fill = "black")),
                   gap = unit(2, "mm")
                 ),
                 top_annotation = upset_top_annotation(m, ylim = c(0,250) ),
))

od = column_order(ht)
cs = comb_size(m)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom",gp = gpar(fontsize = 16))
})
dev.off()


# List of genes shared
sh_lt <- lapply( names(cs[od]), function(i){extract_comb(m, i)})
names(sh_lt) <- c("x1","x2","x3","x4","x5")



