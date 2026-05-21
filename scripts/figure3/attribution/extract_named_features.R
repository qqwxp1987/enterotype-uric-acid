# Extract species attribution data for features explicitly named in Section 3.3
suppressMessages(library(dplyr))

p <- readRDS("results/tables/section_3_7/pathway_contribution_data.rds")
e <- readRDS("results/tables/section_3_7/ec_contribution_data.rds")

cat("=== TYPE 1 PATHWAYS (named in Section 3.3 text) ===\n")
t1_ids <- c("ARGSYN-PWY","ARGSYNBSUB-PWY","GLUTORN-PWY","PWY-7234","PWY-6608")
t1 <- p[["contrib_df"]] %>%
  filter(pathway_id %in% t1_ids) %>%
  group_by(pathway_id, Enterotype, species) %>%
  summarise(mean_pct = mean(contribution_pct), .groups="drop") %>%
  group_by(pathway_id, Enterotype) %>%
  slice_max(mean_pct, n=3) %>%
  filter(Enterotype %in% c(1,4)) %>%
  arrange(pathway_id, Enterotype, desc(mean_pct))
print(as.data.frame(t1), row.names=FALSE)

cat("\n=== TYPE 4 PATHWAYS (named in Section 3.3 text) ===\n")
t4_ids <- c("PWY-5695","PWY-6122","PWY-7219","PEPTIDOGLYCANSYN-PWY","BRANCHED-CHAIN-AA-SYN-PWY")
t4 <- p[["contrib_df"]] %>%
  filter(pathway_id %in% t4_ids) %>%
  group_by(pathway_id, Enterotype, species) %>%
  summarise(mean_pct = mean(contribution_pct), .groups="drop") %>%
  group_by(pathway_id, Enterotype) %>%
  slice_max(mean_pct, n=3) %>%
  filter(Enterotype %in% c(1,4)) %>%
  arrange(pathway_id, Enterotype, desc(mean_pct))
print(as.data.frame(t4), row.names=FALSE)

cat("\n=== DIVERGENT ECs (named in Section 3.3 text) ===\n")
named_ecs <- c("1.4.7.1","3.4.11.5","2.6.1.11","2.4.2.9")
ec_avail <- unique(e[["contrib_long"]]$ec_id)
cat("Available ECs in data:", paste(ec_avail, collapse=", "), "\n")
cat("Requested ECs:", paste(named_ecs, collapse=", "), "\n")
cat("Matched:", paste(intersect(named_ecs, ec_avail), collapse=", "), "\n\n")

ec_named <- e[["contrib_long"]] %>%
  filter(ec_id %in% intersect(named_ecs, ec_avail)) %>%
  group_by(ec_id, Enterotype, species) %>%
  summarise(mean_pct = mean(contribution_pct), .groups="drop") %>%
  group_by(ec_id, Enterotype) %>%
  slice_max(mean_pct, n=3) %>%
  filter(Enterotype %in% c(1,4)) %>%
  arrange(ec_id, Enterotype, desc(mean_pct))
print(as.data.frame(ec_named), row.names=FALSE)

# Also show ALL 8 target ECs from the previous analysis
cat("\n=== ALL 8 TARGET ECs (from Step 3) ===\n")
all_target <- e[["contrib_long"]] %>%
  group_by(ec_id, Enterotype, species) %>%
  summarise(mean_pct = mean(contribution_pct), .groups="drop") %>%
  group_by(ec_id, Enterotype) %>%
  slice_max(mean_pct, n=1) %>%
  filter(Enterotype %in% c(1,4)) %>%
  arrange(ec_id, Enterotype, desc(mean_pct))
print(as.data.frame(all_target), row.names=FALSE)
