/* Documentation for {langev.c} */
/* Last edited on 2008-01-10 10:30:24 by stolfi */ 

#ifndef langev_info_H
#define langev_info_H

#define PROG_HELP \
  PROG_NAME "  \\\n" \
  "  -generations {NGENS} [ -step {GEN_STEP} ] \\\n" \
  "  [ -size {NX} {NY} ] [ -relief {ALTFILE} {LEVEL} ] \\\n" \
  "  -start {IX} {IY} \\\n" \
  "  -prefix {PREFIX}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC_MAIN "\n" \
  "\n" \
  "DEMOGRAPHIC SIMULATION\n" \
  PROG_INFO_DESC_DEMO "\n" \
  "\n" \
  "GEOGRAPHIC DIFFUSION\n" \
  PROG_INFO_DESC_MOVE "\n" \
  "\n" \
  "LANGUAGE EVOLUTION\n" \
  PROG_INFO_DESC_LANG "\n" \
  "\n" \
  "GENOME EVOLUTION\n" \
  PROG_INFO_DESC_GENE "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS"\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnm(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2006-03-15 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Revised on 2006-06-15 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " langev_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC_MAIN \
  "  This program simulates the evolution of natural languages, " \
  " according to a simple model that incorporates many phenomena" \
  " conjectured to occur in the real world.  It also simulates" \
  " the genetic evolution of the population.\n" \
  "\n" \
  "Simulated world\n" \
  "\n" \
  langev_geog_world_INFO "\n" \
  "\n" \
  "Coordinate system\n" \
  "\n" \
  langev_geog_coords_INFO "\n" \
  "\n" \
  "Relief\n" \
  "\n" \
  langev_geog_relief_INFO "\n" \
  "\n" \
  "World state\n" \
  "\n" \
  langev_base_frames_INFO "\n" \
  "\n" \
  "Generation procedure and epochs\n" \
  "\n" \
  "  The simulation consists of a /generation procedure/ that is" \
  " repeated periodically at certain specific times, the /epochs/.  The" \
  " interval between epochs is supposedly the average time between human" \
  " generations.  Each run of the generation procedure takes an" \
  " /old state/ {L,G} and produces a /new state/ {L\',G\'}.\n" \
  "\n" \
  "  The {G} and {L} values evolve by migration, inheritance," \
  " mixing, and random mutation, in ways that will be detailed" \
  " later on.  The {G} values do not influence the simulation, but" \
  " the {L} values (which are taken as surrogates of peple's" \
  " culture and nationality) has a pervasive influence.\n" \
  "\n" \
  "Program output\n" \
  "\n" \
  "  The generation procedure is performed by a user-specified" \
  " number {NGENS}, which are written out as separate PPM" \
  " files.  In these files, water regions" \
  " are displayed in shades of blue, and deserted dry land is painted in" \
  " shades of grey (lighter = higher).  Populated regions are shown in" \
  " various other shades of color, as described" \
  " in the COLOR ENCODING section below.\n" \
  "\n" \
  "  The output file names are \"{PREFIX}-{NNNNNN}-L.ppm}\" and" \
  " \"{PREFIX}-{NNNNNN}-G.ppm}\" where {NNNNNN} is the" \
  " epoch index, from 0 to {NGENS}. If the" \
  " \"-step {GEN_STEP}\" option is given," \
  " the outputimages are produced only for those" \
  " epochs whose indices are multiple of {GEN_STEP}."

#define PROG_INFO_DESC_DEMO \
  "\n" \
  "Population\n" \
  "\n" \
  langev_base_family_INFO "\n" \
  "\n" \
  "Generation procedure\n" \
  "\n" \
  "  At each epoch, when the generation procedure starts, all inhabitants" \
  " of the world are assumed to have just reached adulthood.   Men and women from" \
  " different sites are then /married/ in pairs.  Each marriage produces a couple" \
  " of twin children, one male and one female.  These children learn the" \
  " local language, with some variations, and eventually settle down in" \
  " a nearby site.   These children constitute the" \
  " new state {L\',G\'} of the world.\n" \
  "\n" \
  "Travel\n" \
  "\n" \
  langev_move_travel_INFO "\n" \
  "\n" \
  langev_move_stdiff_INFO "\n" \
  "\n" \
  "Marriages\n" \
  "\n" \
  "  To simulate the marriages, each occupied site is considered in turn, in" \
  " random order. The man living in that site is randomly married to" \
  " one or more women in nearby sites.  (Actually the same model can also be" \
  " interpreted in terms of home-bound men sought by mobile women, if it" \
  " so pleases the reader.)\n" \
  "\n" \
  "  The marriages of a man are restricted to women on sites that can be" \
  " reached from the man's home.  Each wife is drawn from this set, with" \
  " a probability distribution {S} that is biased towards women whose" \
  " home is more easily reached from the man's home, and who speak" \
  " languages closer to his own.  For this purpose, travel over" \
  " shallow water is considered rather difficult.\n" \
  "\n" \
  "  The man's sister is also included in" \
  " the set {S}, with low weight, in order to prevent the set from being" \
  " empty.\n" \
  "\n" \
  "Birth and death\n" \
  "\n" \
  "  The twin children that result from each marriage" \
  " are nominally born on a dry land site which is the" \
  " home of one of the parents.  They will travel and" \
  " settle together, until the next epoch.\n" \
  "\n" \
  "  The children may die before being able to settle down," \
  " (see the Migration section below).  If they succeed in" \
  " settling down, they will" \
  " survive until the next epoch.  By then, all people from" \
  " the previous generation are assumed to have died."

#define PROG_INFO_DESC_MOVE \
  "\n" \
  "Major migration\n" \
  "\n" \
  "  After learning their language, each pair of twins migrates from" \
  " their birth place {pb} to a new dry land site {pc}.  The first stage" \
  " of the migration is in the general direction of a /displacement" \
  " vector/ {d}, to be explained later.\n" \
  "\n" \
  "  In this first stage of the migration, the children will" \
  " travel only on dry land or shallow water, until the total remoteness" \
  " of the path is roughly equal to the Euclidean length of {d}.  For" \
  " this purpose, travel over shallow water is assumed to be quite" \
  " a bit easier than over dry land.  However, the final site {pc}" \
  " will always be on dry land.\n" \
  "\n" \
  "Settling down\n" \
  "\n" \
  "  If the tentative settling position {pd} is already occupied by" \
  " a couple of the same generation, the new couple" \
  " either displaces the previously settled occupants, or /bounces/ leaving those" \
  " occupants undisturbed.  In either case, the couple which is still" \
  " homeless (either because they bounced or because they were" \
  " displaced) gets tentatively settled at some nearby position {pd\'}, reachable" \
  " from {pd} by dry land or shallow water.  The process then repeats: if {pd\'} is" \
  " occupied the homeless couple either displaces the occupants or" \
  " bounces again.  The process ends when a vacant site is" \
  " found, where the homeless couple can settle down; or" \
  " after a fixed number of trials, at which point the unlucky" \
  " children are eaten by wild wolves.\n" \
  "\n" \
  "  To compensate for the infant mortality in the settling process, the" \
  " average fertility over all peoples and times must be kept somewhat" \
  " greater than 1.\n" \
  "\n" \
  "Abduction\n" \
  "\n" \
  "  Rarely, a couple of young twins who have just left the paternal home" \
  " are abducted by the Roc bird and dropped in an entirely random spot" \
  " {pd} of the planet. If {pd} is on water, the children drown; otherwise" \
  " they try to settle there, by the process explained above."

#define PROG_INFO_DESC_LANG \
  "\n" \
  "Learning the language\n" \
  "\n" \
  "  The evolution of languages is simulated as follows.  For each" \
  " newborn couple, a /tutoring family/ is" \
  " chosen at random among the families living on sites reachable from" \
  " their nominal birth site {pb}. The probability" \
  " distribution {S\'} is biased towards" \
  " families whose home is more easily reached from {pb}, and whose" \
  " language is most sililar to that spoken at {pb}.\n" \
  "\n" \
  "  By the time the children are grown up, their language will" \
  " be mostly that of the tutoring family, with a small admixture" \
  " of their mother\'s language, plus a small amount of random" \
  " mutations.\n" \
  "\n" \
  "  The justification for this process is the observation that children" \
  " in the real world learn their ``mother language'' from all their social" \
  " contacts, not just from their parents;" \
  " and that contacts are more likely to occur with people who" \
  " speak languages similar to those of their parents.\n" \
  "\n" \
  "Migration trends\n" \
  "\n" \
  "  The spatial distribution of languages in the real" \
  " world has been significantly affected by the more or less coordinated" \
  " migration of tribes and nations.  The program tries" \
  " to simulate this phenomenon, by assuming that people who" \
  " speak similar languages are more likely to belong to the" \
  " same nation, and therefore are likely to migrate towards similar directions" \
  " and with similar speeds.\n" \
  "\n" \
  "  Specifically, the displacement vector {d} that influences the migration" \
  " of newborn children is defined as a /trend/ vector {t} plus a" \
  " small random vector {r}. The trend vector depends strongly on the" \
  " language {v} spoken by the twins, but" \
  " varies randomly but slowly over time; so that all" \
  " persons speaking similar languages will tend to drift" \
  " together over several successive generations.\n" \
  "\n" \
  "Population prestige\n" \
  "\n" \
  "  Another important factor in the evolution of" \
  " real world languages are the cultural, technlogical," \
  " economic, or military differences between nations," \
  " that often determine the fertility/mortality rate" \
  " of their inhabitants, and which language will" \
  " prevail in bilingual contexts.\n" \
  "\n" \
  "   To simulate these factors, the program defines for each epoch {t}" \
  " and each language {v} a quantity {prs(t,v)}, here called the" \
  " language's /prestige/.  During the generation procedure, the number" \
  " of marriages simulated for each man in epoch {t} is a random" \
  " quantity whose expected value is a motononic funcion of {prs(t,v)}," \
  " where {v} is the man's language.  As in" \
  " the case of migration, the program assumes that speakers of" \
  " similar languages are more likely to belong to the same nation, and" \
  " therefore to have similar nation-related handicaps.   Like the trend vector, the" \
  " prestige is a highly variable function of the language, which" \
  " varies randomly but slowly over time.\n" \
  "\n" \
  "  The prestige of a language is also relevant to the process of" \
  " settling a new couple of twins in the new world state.  When there" \
  " is a clash with a previously settled couple, the probability of" \
  " the outcome (bounce or displace) depend on the prestiges" \
  " of their languages in the new epoch."

#define PROG_INFO_DESC_GENE \
  "\n" \
  "Genetic evolution\n" \
  "\n" \
  "  The genetic makeup of the twin children born of a marriage" \
  " is represented by a binary string, divided into three parts.  One" \
  " part {Gm} (representing the mitochondrial DNA) is inherited from" \
  " the mother only.  The other part {Gf} (representing the male Y" \
  " chromosome) is inherited from the father only).  The bits of the" \
  " third part {Gx} are randomly inherited from either parent," \
  " independently, with equal probabilities.\n" \
  "\n" \
  "  The current version of the" \
  " program does not try to model the diploid character of the" \
  " human genome.  In particular, it assumes that the two X" \
  " chromosomes of the woman are pre-mixed so only" \
  " one copy needs to be recorded.\n" \
  "\n" \
  "  Unlike the language, the genetic part of the state does" \
  " not have any influence on the simulation process, except" \
  " in the genetics of the next generation.  This seems to be" \
  " in accord with the real world: in the long term, genetic" \
  " differences are much less important for linguistic" \
  " and demographic evolution than national/political/cultural factors."

#define PROG_INFO_OPTS \
  "  -generations {NGENS}\n" \
  "    This mandatory directive specifies the number of generations" \
  " to simulate.  The number {NGENS} must be non-nengative.\n" \
  "\n" \
  "  -step {GEN_STEP}\n" \
  "    This optional directive specifies that output images" \
  " should be produced only for epochs whose index is a multiple" \
  " of {GEN_STEP}.  The default is 1 (output all epochs).\n" \
  "\n" \
  "  -relief {ALTFILE} {LEVEL}\n" \
  "    This directive specifies that the dimensions and physical" \
  " geography of the world are given by the file \"{ALTFILE}\"," \
  " which must contain a monochrome image in PGM format.  The" \
  " {LEVEL} parameter is the pixel value in {ALTFILE} that" \
  " corresponds to shallow water.  The user must" \
  " specify either \"-relief\" or \"-size\".\n" \
  "\n" \
  "  -size {NX} {NY}\n" \
  "    This directive specifies that the entire world is a" \
  " featureless dry-land plain, with {NX} columns and {NY}" \
  " rows of sites, all at the same altitude.  The user must" \
  " specify either \"-size\" or \"-relief\"; if both are given," \
  " the relief image {ALTFILE} must have {NX} columns and {NY} rows.\n" \
  "\n" \
  "  -start {IX} {IY}\n" \
  "    This mandatory directive specifies the location of the" \
  " initial couple in epoch 0. The indices range from 0 to {NX-1}" \
  " and from 0 to {NY-1}, respectively, where {NX} and {NY} are" \
  " the counts of site columns and row in the world. \n" \
  "\n" \
  "  -prefix {PREFIX}\n" \
  "    This mandatory directive specifies the common prefix" \
  " for the names of all output files."

#endif
