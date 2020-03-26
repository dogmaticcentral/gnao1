
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
   `drugbank_id`  varchar(10) NOT NULL,
   `pubchem`   int,
   `synonyms`  text CHARACTER SET utf8mb4 ,
   `products`  text CHARACTER SET utf8mb4,
   `targets`   text CHARACTER SET utf8mb4,
   `brands`    text CHARACTER SET utf8mb4,
   `is_prodrug_of` text CHARACTER SET utf8mb4,
    PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


-- to get around importing problem with mysql
-- sudo mv bindingdb_ki.tsv /var/lib/mysql-files/
-- sudo  mysqlimport gnao1  /var/lib/mysql-files/bindingdb_ki.tsv
DROP TABLE IF EXISTS `bindingdb_ki`;
CREATE TABLE `bindingdb_ki`(
   `id` int  NOT NULL  AUTO_INCREMENT,
   `drug_name` text CHARACTER SET utf8mb4,
   `target_symbol`  varchar(50) NOT NULL,
   `ki_nM`  int,
   `mode_of_action` varchar(50) ,
     PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS `guidetopharm_ki`;
CREATE TABLE `guidetopharm_ki`(
   `id` int  NOT NULL  AUTO_INCREMENT,
   `drug_name` text CHARACTER SET utf8mb4,
   `target_symbol`  varchar(100) NOT NULL,
   `ki_nM`  int,
   `mode_of_action` varchar(50) ,
     PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `pdsp_ki`;
CREATE TABLE `pdsp_ki`(
   `id` int  NOT NULL,
   `drug_name` text CHARACTER SET utf8mb4,
   `target_symbol`  varchar(50) NOT NULL,
   `ki_nM`  int,
   `mode_of_action` varchar(50) ,
   `species`  varchar(50) NOT NULL,

     PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `pubchem_ki`;
CREATE TABLE `pubchem_ki`(
   `id` int  NOT NULL  AUTO_INCREMENT,
   `drug_name` text CHARACTER SET utf8mb4,
   `target_symbol`  varchar(50) NOT NULL,
   `ki_nM`  int,
   `mode_of_action` varchar(50) ,
     PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



DROP TABLE IF EXISTS `literature_ki`;
CREATE TABLE `literature_ki`(
   `id` int  NOT NULL  AUTO_INCREMENT,
   `drug_name` text CHARACTER SET utf8mb4,
   `target_symbol`  varchar(50) NOT NULL,
   `ki_nM`  int,
   `mode_of_action` varchar(50) ,
   `pubmed`    int,
   `pubmedcentral`    varchar(20),
   `other_sources` text,
    PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


