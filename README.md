# ThesisCode
```{r - formatting clogit}
#CORRELATIONS FOR ALL CONTINUOUS DATA
data <- read.csv("/Users/alyxhutcheon/Library/CloudStorage/OneDrive-Nexus365/DISS/Data/2. R_CODE/4. Data/retrieval/correlations5.csv")

cor_matrix <- cor(data)
print (cor_matrix)
corrplot(cor_matrix, method = "circle")

```{r - formatting clogit}
#CONDITIONAL LOGISTIC REGRESSION FORMATTING
library(survival)
library(ggplot2)
library(ggeffects)
library(corrplot)
library(emmeans)
library(effects)
library(tidyverse)
library (dplyr)
library(lme4)
library(lmerTest)
library (car)
library(viridis)


adata=read.table(file="/Users/alyxhutcheon/Library/CloudStorage/OneDrive-Nexus365/DISS/Data/2. R_CODE/4. Data/retrieval/retrieval7.csv", header=TRUE, sep=",")
head(adata)

attach(adata)
 
#always start with a plot
 
 
head(adata)

#biomass
adata$TotalBiomass <- (Cattle_pres*175)+(Sheep_pres*24)+(Goat_pres*24)+(Donkey_pres*137.5)

adata

# Assuming adata is your data frame
adata <- adata %>%
  mutate(
    WB_Dist = ifelse(Attack.non == 1, WB_Dist / 10, WB_Dist * 10),
    Urban_Dist = ifelse(Attack.non == 0, Urban_Dist * 10, Urban_Dist),
    NP_Dist = ifelse(Attack.non == 0, NP_Dist * 10, NP_Dist)
  )

# Print the modified data frame
print(adata)

summary2<-adata%>%
         group_by(Ref,Boma_ID, Attack.non, Predator_species,Livestock_species, Number_people, Number_Dogs, Bush_Dist, Outer_boma, Boma_Dist, WB_Dist, NP_Dist, Urban_Dist, Material_1, Total_Presence, TotalBiomass, BomaHeight, BomaWidth, BomaStems, BomaWeakness)%>%
         summarise(total_killed=sum(Livestock_killed)) %>%
         filter(Attack.non == 0 | (Attack.non == 1 & total_killed > 0))

summary2_duplicates <- summary2 %>%
  group_by(Ref, Livestock_species) %>%
  filter(n() > 1)


valid_pairs <- summary2 %>%
  group_by(Ref, Livestock_species) %>%
  filter(n_distinct(Attack.non) == 2) %>%
  ungroup()

summary2

#Filter the dataset to include only the valid pairs
filtered_data <- summary2 %>%
  semi_join(valid_pairs, by = c("Ref", "Livestock_species"))

#For each Ref and Livestock_species, keep one attack and one non-attack record
final_filtered <- filtered_data %>%
  group_by(Ref, Livestock_species, Attack.non) %>%
  slice(1) %>%  # Keep only one record per combination of Ref, Livestock_species, and Attack.non
  ungroup()

# Print the final filtered data
print(final_filtered)

attack_records <- final_filtered %>% filter(Attack.non == 1)
non_attack_records <- final_filtered %>% filter(Attack.non == 0)

#Create a mapping from Ref to Predator_species
predator_species_map <- attack_records %>%
  select(Ref, Predator_species) %>%
  distinct()

#Update non-attack records
non_attack_updated <- non_attack_records %>%
  left_join(predator_species_map, by = "Ref", suffix = c("", ".attack")) %>%
  mutate(Predator_species = Predator_species.attack) %>%
  select(-Predator_species.attack)  # Remove the extra column if not needed

#Combine the updated non-attack records with attack records
final_data <- bind_rows(attack_records, non_attack_updated)

# Print the final dataset
print(final_data)

#Order the final dataset by Ref
final_data_ordered <- final_data %>%
  arrange(Ref)%>%

# Print the final ordered dataset
print(final_data_ordered)

final_data_ordered <- final_data_ordered %>%
  mutate(Ref_num = as.numeric(factor(Ref)))

# View the updated dataset
print(final_data_ordered)

```

```{r - scaling}
#CONDITIONAL LOGISTIC REGRESSION SCALING
final_data_ordered$Number_Dogs <- scale(final_data_ordered$Number_Dogs)
final_data_ordered$Total_Presence <- scale(final_data_ordered$Total_Presence)
final_data_ordered$Number_people <- scale(final_data_ordered$Number_people)
final_data_ordered$Bush_Dist <- scale(final_data_ordered$Bush_Dist)
final_data_ordered$Boma_Dist <- scale(final_data_ordered$Boma_Dist)
final_data_ordered$TotalBiomass <- scale(final_data_ordered$TotalBiomass)
final_data_ordered$Urban_Dist <- scale(final_data_ordered$Urban_Dist)
final_data_ordered$WB_Dist <- scale(final_data_ordered$WB_Dist)
final_data_ordered$NP_Dist<- scale(final_data_ordered$NP_Dist)
final_data_ordered$BomaHeight <- scale(final_data_ordered$BomaHeight)
final_data_ordered$BomaWidth <- scale(final_data_ordered$BomaWidth)

final_data_ordered$BomaStems <- factor(final_data_ordered$BomaStems, levels = c("Few", "Most", "All"))
final_data_ordered$BomaWeakness <- factor(final_data_ordered$BomaWeakness, levels = c("Many", "Some", "Few"))
final_data_ordered$Predator_species <- as.factor(final_data_ordered$Predator_species)
final_data_ordered$Livestock_species <- as.factor(final_data_ordered$Livestock_species)

hydata <- final_data_ordered %>% filter(Predator_species == "Hy")
liondata <- final_data_ordered %>% filter(Predator_species == "Lion")

hydata
liondata

final_data_ordered <- final_data_ordered %>%
  filter(!(Predator_species %in% c("Leopard", "Cheetah")))

final_data_ordered

filtered_attack <- final_data_ordered %>% filter(Attack.non == "1")
filtered_noattack <- final_data_ordered %>% filter(Attack.non == "0")



```

```{r - clogit model}
#CONDITIONAL LOGISTIC REGRESSION - HYAENA AND LION MODELS
hymodel <- clogit(
  Attack.non ~  Number_Dogs + BomaHeight + BomaWidth + BomaStems + BomaWeakness + Outer_boma  + Total_Presence + Number_people + Bush_Dist + NP_Dist + WB_Dist +
  strata(Ref_num),
  data = hydata)

summary(hymodel)

#interpretation for clogit 

anova1 <- anova(hymodel, type="III")
print(anova1)

lionmodel <- clogit(
  Attack.non ~  Number_Dogs + BomaHeight + BomaWidth + BomaWeakness + Outer_boma + 
    Livestock_species + Total_Presence + Number_people + Bush_Dist + NP_Dist + WB_Dist +
  strata(Ref_num),
  data = liondata)

summary(lionmodel)

anova2 <- anova(lionmodel, type="III")
print(anova2)

emm1<- emmeans(lionmodel, ~ BomaWeakness)
tukeylionBW <- pairs(emm1, adjust = "tukey")
print(tukeylionBW)


```

```{r - BomaStems alluvial}
#CONDITIONAL LOGISTIC REGRESSION DESCRIPTIVE DATA - LION ALLUVIALS 
#liondata swapped out for hydata and vice versa
combined_data2 <- liondata %>%
  mutate(BomaType = if_else(Attack.non == 1, "Attack", "NonAttack")) %>%
  select(Ref_num, BomaType, BomaStems)

# Expanding each entry based on BomaType
expanded_data2 <- combined_data2 %>%
  distinct() %>%
  arrange(Ref_num, BomaType, BomaStems) %>%
  group_by(Ref_num, BomaType) %>%
  mutate(entry_id = row_number()) %>%
  ungroup()

#Pivot the data to have separate columns for Attack and Non-Attack
separated_data2 <- expanded_data2 %>%
  pivot_wider(
    names_from = BomaType,
    values_from = BomaStems,
    names_prefix = "BomaStems_",
    values_fill = list(BomaStems = NA)
  ) %>%
  arrange(Ref_num, entry_id) %>%
  select(-entry_id)

print(separated_data2)

cleaned_data2 <- separated_data2 %>%
  replace_na(list(BomaStems_Attack = "Few", BomaStems_NonAttack = "Few"))

# Group by combinations of BomaWeakness values and count occurrences
combination_counts2 <- cleaned_data2 %>%
  group_by(BomaStems_Attack, BomaStems_NonAttack) %>%
  summarize(count = n(), .groups = "drop")

combination_counts2 <- combination_counts2 %>%
  mutate(
    BomaStems_Attack = factor(BomaStems_Attack, levels = c("All", "Most", "Few")),
    BomaStems_NonAttack = factor(BomaStems_NonAttack, levels = c("All", "Most", "Few"))
  )

# Inspect the frequency of each combination
print(combination_counts2)

colors <-  viridis(3, begin = 0.2, end = 1)

alluvial_plot2 <- ggplot(combination_counts2,
                        aes(axis1 = BomaStems_Attack,
                            axis2 = BomaStems_NonAttack,
                            y = count)) +
  geom_alluvium(aes(fill = BomaStems_Attack), width = 1/12) +
  geom_stratum(width = 1/8, fill = "grey", color = "black") +
  geom_text(stat = "stratum",
            aes(label = ifelse(after_stat(count) > min_size, as.character(after_stat(stratum)), "")),
            color = "black", size = 5, hjust = 0.5, vjust = 0.4, show.legend = FALSE, fontface = "bold") +
  scale_x_discrete(name = NULL,
                   labels = c("Attack", "Non Attack"),
                   limits = c("BomaStems_Attack", "BomaStems_NonAttack"),
                   expand = c(.05, .05),
                   position = "top", position_dodge(width=0.8)) +
  scale_fill_manual(values=colors) +
  labs(title = "",
       x = "Boma Stems",
       y = "Count",
       fill = "Boma Stems Out") +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", vjust = 2),
    axis.title.x = element_text(size = 14, face = "bold", vjust = -2),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Print the plot
print(alluvial_plot2)

ggsave(alluvial_plot2, filename = "totalalluvial_BS.png", width = 2000, height = 1000, units = "px")


```

```{r - categoric plots}
#CONDITIONAL LOGISTIC REGRESSION DESCRIPTIVE DATA - CATEGORICAL 
# Ensure correct library is loaded
library(dplyr)
library(ggplot2)
library(viridis)

# Calculate the count of occurrences for each combination of "Predator_species" and "Livestock_species"
# Ensure correct library is loaded
library(dplyr)
library(ggplot2)
library(viridis)

# Calculate the count of occurrences for each combination of "Predator_species" and "Livestock_species"
data_counts1 <- filtered_attack %>%
  group_by(Predator_species, Livestock_species) %>%
  summarize(Count = n(), .groups = 'drop')

# Calculate the total count of each "Predator_species"
total_counts1 <- data_counts1 %>%
  group_by(Predator_species) %>%
  summarise(Total = sum(Count), .groups = 'drop')

# Calculate proportions for each "Livestock_species" within each "Predator_species"
data_counts1 <- data_counts1 %>%
  left_join(total_counts1, by = "Predator_species") %>%
  mutate(Proportion = Count / Total)

# Create a new data frame with total sample size for each predator species
total_sample_size1 <- data_counts1 %>%
  group_by(Predator_species) %>%
  summarise(TotalSampleSize = sum(Count, na.rm = TRUE), .groups = 'drop')

# Merge the total sample size information into the original data frame
data_counts1 <- merge(data_counts1, total_sample_size1, by = "Predator_species")

# Update New_Name with the correct format and "n="
data_counts1$New_Name <- paste(data_counts1$Predator_species, "n=", data_counts1$TotalSampleSize, sep = "")

# Replace specific labels with 'Lion' and 'Hyaena' and add new names
data_counts1 <- data_counts1 %>%
  mutate(New_Name = case_when(
    Predator_species == "Lion" ~ paste("Lion\nn=", TotalSampleSize),
    Predator_species == "Hyaena" ~ paste("Hyaena\nn=", TotalSampleSize),
    TRUE ~ New_Name
  ))

# Check the data to confirm the labels
print(unique(data_counts1$Predator_species))

# Create a threshold for displaying text labels (e.g., only show text on sections with proportion > 0.1)
threshold <- 0.1

# Define color palette
colour_palette <- viridis(3, begin = 0.2, end = 1)

# Create the plot
categories <- ggplot(data_counts1, aes(x = New_Name, y = Proportion, fill = Livestock_species)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ifelse(Proportion > threshold, paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  scale_fill_manual(values = colour_palette,
                    labels = c("Cattle", "Shoat")) +
  labs(x = "", y = "Proportion (%)", fill = "Livestock species") +
  theme(text = element_text(family="Calibri"),
        axis.text.y = element_text(size = 2)) + 
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 13))+
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) +  
  theme(legend.text = element_text(face = "bold", size = 10))

# Display the plot
print(categories)

# Save the plot if needed

# Save the plot if needed
ggsave(categories, filename = "LivestockxPred.png", width = 2000, height = 1000, units = "px")

#Calculate the count of occurrences for each combination of "Predatorspecies" and "BomaWeakness species"
data_counts2 <- filtered_attack %>%
  group_by(Predator_species, BomaWeakness) %>%
  summarize(Count = n())

# Calculate the total count of each "DataSource"
total_counts2 <- data_counts2 %>%
  group_by(Predator_species) %>%
  summarise(Total = sum(Count))

# Calculate proportions for each "Category" within each "DataSource"
data_counts2 <- data_counts2 %>%
  left_join(total_counts2, by = "Predator_species") %>%
  mutate(Proportion = Count / Total)

# Create a new data frame with total sample size for each country
total_sample_size2 <- data_counts2 %>%
  group_by(Predator_species) %>%
  summarise(TotalSampleSize = sum(Count, na.rm = TRUE))

# Merge the total sample size information into the original data frame
data_counts2 <- merge(data_counts2, total_sample_size2, by = "Predator_species")
data_counts2$New_Name <- paste (data_counts2$Predator_species, data_counts2$TotalSampleSize, sep = " n=")

# Example assuming data_counts8 is your data frame
data_counts2$New_Name <- gsub(" survey", "", data_counts2$New_Name, ignore.case = TRUE)

# Insert line breaks in country names at a specific point (e.g., after the first space)
data_counts2$New_Name <- gsub(" n", "\nn", data_counts2$New_Name)

# Create a threshold for displaying text labels (e.g., only show text on sections with proportion > 0.05)
threshold <- 0.1

colour_palette <- viridis(3, begin = 0.2, end = 1)

#create plot of proportions flipped of eg. Welsh people who said somewhere in Wales
categories2 <- ggplot(data_counts2, aes(x = New_Name, y = Proportion, fill = BomaWeakness)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ifelse(Proportion > threshold, paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  scale_fill_manual(values = colour_palette,
                    labels = c("Many", "Some", "Few")) +
  labs(x = "", y = "Proportion (%)", fill = "Boma Weakness") +
  theme(text = element_text(family="Calibri")) + 
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 13))+
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.y= element_text(face = "bold", size = 10))+  
  theme(legend.text = element_text(face= "bold", size = 10))

categories2

ggsave(categories2, filename = "BWxPred.png", width = 2000, height = 1000, units = "px")


#Calculate the count of occurrences for each combination of "Predatorspecies" and "BomaStems"
data_counts3 <- filtered_attack %>%
  group_by(Predator_species, BomaStems) %>%
  summarize(Count = n())

# Calculate the total count of each "DataSource"
total_counts3 <- data_counts3 %>%
  group_by(Predator_species) %>%
  summarise(Total = sum(Count))

# Calculate proportions for each "Category" within each "DataSource"
data_counts3 <- data_counts3 %>%
  left_join(total_counts3, by = "Predator_species") %>%
  mutate(Proportion = Count / Total)

# Create a new data frame with total sample size for each country
total_sample_size3 <- data_counts3 %>%
  group_by(Predator_species) %>%
  summarise(TotalSampleSize = sum(Count, na.rm = TRUE))

# Merge the total sample size information into the original data frame
data_counts3 <- merge(data_counts3, total_sample_size3, by = "Predator_species")
data_counts3$New_Name <- paste (data_counts3$Predator_species, data_counts3$TotalSampleSize, sep = " n=")

# Example assuming data_counts8 is your data frame
data_counts3$New_Name <- gsub(" survey", "", data_counts3$New_Name, ignore.case = TRUE)

# Insert line breaks in country names at a specific point (e.g., after the first space)
data_counts3$New_Name <- gsub(" n", "\nn", data_counts3$New_Name)

# Create a threshold for displaying text labels (e.g., only show text on sections with proportion > 0.05)
threshold <- 0.1

colour_palette <- viridis(3, begin = 0.2, end = 1)

#create plot of proportions flipped of eg. Welsh people who said somewhere in Wales
categories3 <- ggplot(data_counts3, aes(x = New_Name, y = Proportion, fill = BomaStems)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ifelse(Proportion > threshold, paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  scale_fill_manual(values = colour_palette,
                    labels = c("Few", "Most", "All")) +
  labs(x = "", y = "Proportion (%)", fill = "Boma Stems") +
  theme(text = element_text(family="Calibri")) + 
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 13))+
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.y= element_text(face = "bold", size = 10))+  
  theme(legend.text = element_text(face= "bold", size = 10))

categories3

ggsave(categories3, filename = "BomaStemsxPred.png", width = 2000, height = 1000, units = "px")

#Calculate the count of occurrences for each combination of "Predatorspecies" and "Livestock species"
data_counts4 <- filtered_attack %>%
  group_by(Attack.non, Livestock_species) %>%
  summarize(Count = n())

# Calculate the total count of each "DataSource"
total_counts4 <- data_counts4 %>%
  group_by(Attack.non) %>%
  summarise(Total = sum(Count))

# Calculate proportions for each "Category" within each "DataSource"
data_counts4 <- data_counts4 %>%
  left_join(total_counts4, by = "Attack.non") %>%
  mutate(Proportion = Count / Total)

# Create a new data frame with total sample size for each country
total_sample_size4 <- data_counts4 %>%
  group_by(Attack.non) %>%
  summarise(TotalSampleSize = sum(Count, na.rm = TRUE))

# Merge the total sample size information into the original data frame
data_counts4 <- merge(data_counts4, total_sample_size4, by = "Attack.non")
data_counts4$New_Name <- paste (data_counts4$Attack.non, data_counts4$TotalSampleSize, sep = " n=")

# Example assuming data_counts8 is your data frame
data_counts4$New_Name <- gsub(" survey", "", data_counts4$New_Name, ignore.case = TRUE)

# Insert line breaks in country names at a specific point (e.g., after the first space)
data_counts4$New_Name <- gsub(" n", "\nn", data_counts4$New_Name)

# Create a threshold for displaying text labels (e.g., only show text on sections with proportion > 0.05)
threshold <- 0.1

colour_palette <- viridis(3, begin = 0.2, end = 1)

#create plot of proportions flipped of eg. Welsh people who said somewhere in Wales
categories4 <- ggplot(data_counts4, aes(x = New_Name, y = Proportion, fill = Livestock_species)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ifelse(Proportion > threshold, paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  scale_fill_manual(values = colour_palette,
                    labels = c("Cattle", "Shoat")) +
  labs(x = "", y = "Proportion (%)", fill = "Livestock species") +
  theme(text = element_text(family="Calibri")) + 
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 13))+
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.y= element_text(face = "bold", size = 10))+  
  theme(legend.text = element_text(face= "bold", size = 10))

categories4

ggsave(categories4, filename = "LivestockxAttacknon.png", width = 2000, height = 1000, units = "px")

#Calculate the count of occurrences for each combination of "Predatorspecies" and "Livestock species"
data_counts5 <- filtered_attack %>%
  group_by(Attack.non, BomaWeakness) %>%
  summarize(Count = n())

# Calculate the total count of each "DataSource"
total_counts5 <- data_counts5 %>%
  group_by(Attack.non) %>%
  summarise(Total = sum(Count))

# Calculate proportions for each "Category" within each "DataSource"
data_counts5 <- data_counts5 %>%
  left_join(total_counts5, by = "Attack.non") %>%
  mutate(Proportion = Count / Total)

# Create a new data frame with total sample size for each country
total_sample_size5 <- data_counts5 %>%
  group_by(Attack.non) %>%
  summarise(TotalSampleSize = sum(Count, na.rm = TRUE))

# Merge the total sample size information into the original data frame
data_counts5 <- merge(data_counts5, total_sample_size5, by = "Attack.non")
data_counts5$New_Name <- paste (data_counts5$Attack.non, data_counts5$TotalSampleSize, sep = " n=")

# Example assuming data_counts8 is your data frame
data_counts5$New_Name <- gsub(" survey", "", data_counts5$New_Name, ignore.case = TRUE)

# Insert line breaks in country names at a specific point (e.g., after the first space)
data_counts5$New_Name <- gsub(" n", "\nn", data_counts5$New_Name)

# Create a threshold for displaying text labels (e.g., only show text on sections with proportion > 0.05)
threshold <- 0.1

colour_palette <- viridis(3, begin = 0.2, end = 1)

#create plot of proportions flipped of eg. Welsh people who said somewhere in Wales
categories5 <- ggplot(data_counts5, aes(x = New_Name, y = Proportion, fill = BomaWeakness)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ifelse(Proportion > threshold, paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  scale_fill_manual(values = colour_palette,
                    labels = c("Many", "Some", "Few")) +
  labs(x = "", y = "Proportion (%)", fill = "Boma Weakness") +
  theme(text = element_text(family="Calibri")) + 
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 13))+
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.y= element_text(face = "bold", size = 10))+  
  theme(legend.text = element_text(face= "bold", size = 10))

categories5

ggsave(categories5, filename = "BWxAttacknon.png", width = 2000, height = 1000, units = "px")


#Calculate the count of occurrences for each combination of "Predatorspecies" and "Livestock species"
data_counts6 <- filtered_attack %>%
  group_by(Attack.non, BomaStems) %>%
  summarize(Count = n())

# Calculate the total count of each "DataSource"
total_counts6 <- data_counts6 %>%
  group_by(Attack.non) %>%
  summarise(Total = sum(Count))

# Calculate proportions for each "Category" within each "DataSource"
data_counts6 <- data_counts6 %>%
  left_join(total_counts6, by = "Attack.non") %>%
  mutate(Proportion = Count / Total)

# Create a new data frame with total sample size for each country
total_sample_size6 <- data_counts6 %>%
  group_by(Attack.non) %>%
  summarise(TotalSampleSize = sum(Count, na.rm = TRUE))

# Merge the total sample size information into the original data frame
data_counts6 <- merge(data_counts6, total_sample_size6, by = "Attack.non")
data_counts6$New_Name <- paste (data_counts6$Attack.non, data_counts6$TotalSampleSize, sep = " n=")

# Example assuming data_counts8 is your data frame
data_counts6$New_Name <- gsub(" survey", "", data_counts6$New_Name, ignore.case = TRUE)

# Insert line breaks in country names at a specific point (e.g., after the first space)
data_counts6$New_Name <- gsub(" n", "\nn", data_counts6$New_Name)

# Create a threshold for displaying text labels (e.g., only show text on sections with proportion > 0.05)
threshold <- 0.1

colour_palette <- viridis(3, begin = 0.2, end = 1)

#create plot of proportions flipped of eg. Welsh people who said somewhere in Wales
categories6 <- ggplot(data_counts6, aes(x = New_Name, y = Proportion, fill = BomaStems)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = ifelse(Proportion > threshold, paste0(round(Proportion * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3) +
  scale_fill_manual(values = colour_palette,
                    labels = c("Few", "Most", "All")) +
  labs(x = "", y = "Proportion (%)", fill = "Boma Stems") +
  theme(text = element_text(family="Calibri")) + 
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face = "bold", size = 13))+
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.y= element_text(face = "bold", size = 10))+  
  theme(legend.text = element_text(face= "bold", size = 10))

categories6

ggsave(categories6, filename = "StemsxAttacknon.png", width = 2000, height = 1000, units = "px")



```

```{r descriptive liondata}
#CONDITIONAL LOGISTIC REGRESSION DESCRIPTIVE DATA - LION HISTOGRAMS 
#Calculate the count of occurrences for each combination of "Predatorspecies" and "Livestock species"

colour_palette <- viridis(3, begin = 0.2, end = 1)

histogram_plot1 <- ggplot(liondata, aes(x = BomaHeight, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Boma Height (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

# Display the histogram
histogram_plot1

ggsave(histogram_plot1, filename = "LionBH.png", width = 2000, height = 1000, units = "px")


histogram_plot2 <- ggplot(liondata, aes(x = BomaWidth, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Boma Width (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )
histogram_plot2

ggsave(histogram_plot2, filename = "LionBW.png", width = 2000, height = 1000, units = "px")

histogram_plot3 <- ggplot(liondata, aes(x = Total_Presence, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Number of Livestock Present", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )
histogram_plot3

ggsave(histogram_plot3, filename = "LionTP.png", width = 2000, height = 1000, units = "px")

histogram_plot4 <- ggplot(liondata, aes(x = Number_people, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Number_people", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot4

ggsave(histogram_plot4, filename = "LionNPeople.png", width = 2000, height = 1000, units = "px")

histogram_plot5 <- ggplot(liondata, aes(x = Bush_Dist, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Bush Proximity (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot5

ggsave(histogram_plot5, filename = "LionBush.png", width = 2000, height = 1000, units = "px")

histogram_plot6 <- ggplot(liondata, aes(x = NP_Dist, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "NP Proximity (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot6

ggsave(histogram_plot6, filename = "LionNP.png", width = 2000, height = 1000, units = "px")


histogram_plot7 <- ggplot(liondata, aes(x = WB_Dist, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "WB Proximity (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot7

ggsave(histogram_plot7, filename = "LionWB.png", width = 2000, height = 1000, units = "px")

histogram_plot8 <- ggplot(liondata, aes(x = Number_Dogs, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Number of Dogs", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot8

ggsave(histogram_plot8, filename = "LionDogs.png", width = 2000, height = 1000, units = "px")


```

```{r - descriptive hydata}
#CONDITIONAL LOGISTIC REGRESSION DESCRIPTIVE DATA - HYAENA HISTOGRAMS
#Calculate the count of occurrences for each combination of "Predatorspecies" and "Livestock species"

histogram_plot9 <- ggplot(hydata, aes(x = BomaHeight, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Boma Height (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

# Display the histogram
histogram_plot9

ggsave(histogram_plot9, filename = "HyBH.png", width = 2000, height = 1000, units = "px")

histogram_plot10 <- ggplot(hydata, aes(x = BomaWidth, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Boma Width (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )
histogram_plot10

ggsave(histogram_plot10, filename = "HyBW.png", width = 2000, height = 1000, units = "px")


histogram_plot11 <- ggplot(hydata, aes(x = Total_Presence, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Number of Livestock Present", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )
histogram_plot11

ggsave(histogram_plot11, filename = "HyTP.png", width = 2000, height = 1000, units = "px")


histogram_plot12 <- ggplot(hydata, aes(x = Number_people, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Number_people", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot12

ggsave(histogram_plot12, filename = "HyNPeople.png", width = 2000, height = 1000, units = "px")

histogram_plot13 <- ggplot(hydata, aes(x = Bush_Dist, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Bush Proximity (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot13

ggsave(histogram_plot13, filename = "HyBush.png", width = 2000, height = 1000, units = "px")


histogram_plot14 <- ggplot(hydata, aes(x = NP_Dist, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "NP Proximity (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot14

ggsave(histogram_plot14, filename = "HyNP.png", width = 2000, height = 1000, units = "px")


histogram_plot15 <- ggplot(hydata, aes(x = WB_Dist, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "WB Proximity (m)", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot15

ggsave(histogram_plot15, filename = "HyWB.png", width = 2000, height = 1000, units = "px")


histogram_plot16 <- ggplot(hydata, aes(x = Number_Dogs, fill = factor(Attack.non))) +
  geom_histogram(bins = 10, color = "black", position = "dodge") +
  scale_fill_manual(values = colour_palette,   # Define your custom colors here
                    labels = c("Non-attack", "Attack")) +
  labs(x = "Number of Dogs", y = "Count", fill = "Status") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

histogram_plot16

ggsave(histogram_plot16, filename = "HyDogs.png", width = 2000, height = 1000, units = "px")


```

```{r - severity model}
#LOGISTIC MIXED EFFECTS MODEL - SEVERITY MODEL
bdata=read.table(file="/Users/alyxhutcheon/Library/CloudStorage/OneDrive-Nexus365/DISS/Data/2. R_CODE/4. Data/retrieval/retrieval7.csv", header=TRUE, sep=",")
head(bdata)


#biomass
bdata$TotalBiomass <- (bdata$Cattle_pres*175)+(bdata$Shoat_pres*24)+(bdata$Donkey_pres*137.5)

bdata


summary2<-bdata%>%
         group_by(Boma_ID,Attack.non,Predator_species,Livestock_species, Number_people, Number_Dogs, Bush_Dist, Outer_boma, Boma_Dist, WB_Dist, NP_Dist, Urban_Dist, Total_Presence, TotalBiomass, BomaHeight, BomaWidth, BomaStems, BomaWeakness)%>%
         summarise(total_killed=sum(Livestock_killed)) %>%
        filter(total_killed>0)
 
summary2

#severity calc
calculate_severity <- function(Livestock_species, total_killed) {
  multiplier <- switch(Livestock_species,
                       "Cattle" = 175,
                       "Calf" = 175,
                       "Shoat" = 24,
                       "Donkey" = 137.5,
                       1)  # Default multiplier if species not found
  return(total_killed * multiplier)
}
summary3 <- summary2 %>%
  mutate(severity = mapply(calculate_severity, Livestock_species, total_killed))
summary3$severityprop <- summary3$severity/summary2$TotalBiomass
summary3 <- summary3 %>% filter(is.finite(severityprop))

print(summary3)

summary3$Number_Dogs <- scale(summary3$Number_Dogs)
summary3$Total_Presence <- scale(summary3$Total_Presence)
summary3$Number_people <- scale(summary3$Number_people)
summary3$Bush_Dist <- scale(summary3$Bush_Dist)
summary3$Urban_Dist <- scale(summary3$Urban_Dist)
summary3$WB_Dist <- scale(summary3$WB_Dist)
summary3$NP_Dist<- scale(summary3$NP_Dist)
summary3$BomaHeight <- scale(summary3$BomaHeight)
summary3$BomaWidth <- scale(summary3$BomaWidth)

summary3$BomaStems <- factor(summary3$BomaStems, levels = c("Few", "Most", "All"))
summary3$BomaWeakness <- factor(summary3$BomaWeakness, levels = c("Many", "Some", "Few"))
summary3$Predator_species <- as.factor(summary3$Predator_species)
summary3$Livestock_species <- as.factor(summary3$Livestock_species)

hist(summary3$severityprop)

summary3

boma_sev3 <- lmer(severityprop ~ Predator_species * (Number_Dogs+BomaHeight+BomaWidth+BomaStems+BomaWeakness + Outer_boma) + Livestock_species + Total_Presence + Number_people + Bush_Dist + NP_Dist + WB_Dist + (1|Boma_ID), data = summary3)
summary(boma_sev3)


anovamodel<- anova(boma_sev3, type='III')
summary(anovamodel)
print(anovamodel)

emm1<- emmeans(boma_sev3, ~ Predator_species*BomaWeakness)
PSxBW <- pairs(emm1, adjust = "tukey")
print(PSxBW)

emm2 <-emmeans(boma_sev3, ~Predator_species)
PS <- pairs(emm2, adjust ="tukey")
print(PS)

emm3 <- emmeans(boma_sev3, ~ BomaWeakness)
BW <- pairs(emm3, adjust ="tukey")
print(BW)

emm4 <- emmeans(boma_sev3, ~ Livestock_species)
LS <- pairs(emm4, adjust ="tukey")
print(LS)



BW_pred <- ggpredict(boma_sev3, terms = c("BomaWeakness", "Predator_species"), type = "fe", ci.level = 0.95) 

custom_colors <- c("Lion" = "lightsalmon2", "Hy" = "lightskyblue4")
custom_labels <- c("Lion" = "Lion", "Hy" = "Hyaena")

BW_plot <- ggplot(BW_pred)+   
   scale_y_continuous(labels = scales::number_format(), limits = c(0.01, 0.07)) +
  theme_bw()+   
  geom_pointrange(aes(x = x, 
                      y = predicted,  
                      ymin = conf.low,
                      ymax = conf.high, 
                      color = group),  # Use color instead of shape
                  position = position_dodge(width = 0.5)) + 
  labs(title = "", x = "Boma Weakness", y = "Severity") +
  coord_flip() +
  geom_hline(yintercept = 3, colour = gray(1/2), lty = 2) +
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Change to color
  theme(axis.title.x = element_text(face = "bold", size = 13),  
        axis.title.y = element_text(face = "bold", size = 13),  
        axis.text.x = element_text(face = "bold", size = 10),  
        axis.text.y = element_text(face = "bold", size = 10),  
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "right",  # Set legend position
        legend.text = element_text(face = "bold", size = 10),  # Customize legend text
        legend.title = element_text(face = "bold", size = 13)) +  # Customize legend title
  guides(color = guide_legend(title = NULL))  # Customize color guide

# Display the plot
BW_plot

ggsave(BW_plot, filename = "BWPlot.png", width = 1500, height = 1000, units = "px")


LS_pred <- ggpredict(boma_sev3, terms = c("Livestock_species"), type = "fe", ci.level = 0.95)

LS_plot <- ggplot(LS_pred)+   
   scale_y_continuous(labels = scales::number_format(), limits = c(0.001, 0.048)) +
  theme_bw()+   
  geom_pointrange(aes(x=x, 
                      y=predicted,  ymin=conf.low,ymax=conf.high, color=group), position=position_dodge(width = 0.5),)+ 
  labs(title = "", x="Livestock species", y="Severity") +
  coord_flip() +
  #scale_shape_manual(values = c("Lion" = 16, "Hyaena" = 4)) +
  geom_hline(yintercept = 3, colour = gray(1/2), lty = 2)+
  theme(axis.title.x = element_text(face = "bold", size = 13))+  
  theme(axis.title.y = element_text(face = "bold", size = 13))+  
  theme(axis.title.x = element_text(face = "bold"))+  
  theme(axis.text.x= element_text(face = "bold", size = 10))+  
  theme(axis.text.y= element_text(face = "bold", size = 10))+  
 #theme(legend.text = element_text(face= "bold", size = 10))+  
 #theme(legend.title = element_text(face= "bold", size = 13))+  
  theme(plot.title = element_text(face= "bold", size = 13)) +
  theme(legend.position = "none") +
  guides(shape = guide_legend(title = NULL))
 
LS_plot

ggsave(LS_plot, filename = "LSPlot.png", width = 1100, height = 1000, units = "px")


# Generate predictions
TS_Pred <- ggpredict(boma_sev3, terms = c("Total_Presence"), type = "fe", ci.level = 0.95)

# Plotting in linear format
TS_Plot <- ggplot(TS_Pred) +
  scale_y_continuous(labels = scales::number_format(), limits = c(-0.024706746, 0.043)) +
  theme_bw() +
geom_line(aes(x = x, y = predicted, color = group), size = 1) +  # Add lines to show trends
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +  # Add confidence intervals as shaded areas
  labs(title = "", x = "Total_Presence", y = "Severity") +
  # scale_shape_manual(values = c("Lion" = 16, "Hyaena" = 4)) + # Remove if not using shapes
  geom_hline(yintercept = 3, colour = gray(1/2), linetype = "dashed") +
  theme(axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom") +  # Adjust to show legend if needed
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))  # Adjust if you show legend

# Display the plot
TS_Plot

ggsave(TS_Plot, filename = "TSPlot.png", width = 1000, height = 1000, units = "px")

current_directory <- getwd()
print(current_directory)

```

```{r - vif}
#VIFs
vif(hymodel)
vif (lionmodel)
vif(boma_sev3)
```

