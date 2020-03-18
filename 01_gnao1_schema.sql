
DROP TABLE IF EXISTS `cases`;
CREATE TABLE `cases` (
  `id` mediumint(9) NOT NULL AUTO_INCREMENT,
  `cdna`      varchar(100),
  `protein`   varchar(50),
  `pubmed`    int,
  `sex`       char(1),
  `phenotype` text,
  `epilepsy`  tinyint,
  `movement`  tinyint,
  `treatment_effective_E` text,
  `treatment_ineff_E` text,
  `treatment_effective_MD` text,
  `treatment_ineff_MD` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `publications`;
CREATE TABLE `publications` (
   `id` mediumint(9) NOT NULL AUTO_INCREMENT,
   `pubmed`    int,
   `reference`  text CHARACTER SET utf8mb4 NOT NULL,
   `pubmedcentral`    varchar(20),
   `other_xref` text,
    PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



DROP TABLE IF EXISTS `drugs`;
CREATE TABLE `drugs` (
   `id` mediumint(9) NOT NULL AUTO_INCREMENT,
   `name` text CHARACTER SET utf8mb4,
   `pubchem`   int,
   `synonyms`  text CHARACTER SET utf8mb4 ,
   `products`  text CHARACTER SET utf8mb4,
   `targets`   text CHARACTER SET utf8mb4,
   `brands`    text CHARACTER SET utf8mb4,
   `is_prodrug_of` text CHARACTER SET utf8mb4,
    PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



DROP TABLE IF EXISTS `affinities`;
CREATE TABLE `affinities` (
   `id` mediumint(9) NOT NULL,
   `drug_name` text  CHARACTER SET utf8mb4,
   `target_symbol`   varchar(20) NOT NULL,
   `ki`  float,
   `kidb_id` mediumint(9) ,
    PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;




