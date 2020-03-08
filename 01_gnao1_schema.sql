
DROP TABLE IF EXISTS `cases`;
CREATE TABLE `cases` (
  `id` mediumint(9) NOT NULL AUTO_INCREMENT,
  `cdna`      varchar(100) NOT NULL,
  `protein`   varchar(50),
  `pubmed`    mediumint(9),
  `epilepsy`  tinyint,
  `movement`  tinyint,
  `treatment_efficient_E` text,
  `treatment_ineff_E` text,
  `treatment_efficient_MD` text,
  `treatment_ineff_MD` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS `publications`;
CREATE TABLE `publications` (
   `id` mediumint(9) NOT NULL AUTO_INCREMENT,
   `pubmed`    mediumint(9),
   `reference`  text,
   `pubmedcentral`    varchar(20),
   `other_xref` text,
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



DROP TABLE IF EXISTS `drugs`;
CREATE TABLE `drugs` (
   `id` mediumint(9) NOT NULL AUTO_INCREMENT,
   `name`  varchar(100),
   `alt_names`  text,
   `class` text,
   `mechanism_main`  text,
   `mechanism_minor`  text,
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


