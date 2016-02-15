/** @file cmdline_enumeration.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.6
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_ENUMERATION_H
#define CMDLINE_ENUMERATION_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "enumeration"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "enumeration"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "@VERSION_MAJOR@.@VERSION_MINOR@.@VERSION_JUNIOR@"
#endif

enum enum_enumeration { enumeration__NULL = -1, enumeration_arg_ntl = 0, enumeration_arg_babai, enumeration_arg_lp, enumeration_arg_ln };
enum enum_dComp { dComp__NULL = -1, dComp_arg_success = 0, dComp_arg_binary };
enum enum_rComp { rComp__NULL = -1, rComp_arg_length = 0, rComp_arg_piece };

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int dimension_arg;	/**< @brief lattice dimension for LWE samples.  */
  char * dimension_orig;	/**< @brief lattice dimension for LWE samples original value given at command line.  */
  const char *dimension_help; /**< @brief lattice dimension for LWE samples help description.  */
  long modulus_arg;	/**< @brief lattice modulus.  */
  char * modulus_orig;	/**< @brief lattice modulus original value given at command line.  */
  const char *modulus_help; /**< @brief lattice modulus help description.  */
  double sigma_arg;	/**< @brief standard deviaton of error distribution.  */
  char * sigma_orig;	/**< @brief standard deviaton of error distribution original value given at command line.  */
  const char *sigma_help; /**< @brief standard deviaton of error distribution help description.  */
  int beta_arg;	/**< @brief blocksize of BKZ reduction.  */
  char * beta_orig;	/**< @brief blocksize of BKZ reduction original value given at command line.  */
  const char *beta_help; /**< @brief blocksize of BKZ reduction help description.  */
  enum enum_enumeration enumeration_arg;	/**< @brief which enumeration algorithm to use (lp = Lindner Peikert, ln = Liu Nguyen) (default='ln').  */
  char * enumeration_orig;	/**< @brief which enumeration algorithm to use (lp = Lindner Peikert, ln = Liu Nguyen) original value given at command line.  */
  const char *enumeration_help; /**< @brief which enumeration algorithm to use (lp = Lindner Peikert, ln = Liu Nguyen) help description.  */
  double babaiBound_arg;	/**< @brief when running LengthPruning, this factor controlls the generation of the R sequecne and especially when to switch to Babais enumeration (default='4').  */
  char * babaiBound_orig;	/**< @brief when running LengthPruning, this factor controlls the generation of the R sequecne and especially when to switch to Babais enumeration original value given at command line.  */
  const char *babaiBound_help; /**< @brief when running LengthPruning, this factor controlls the generation of the R sequecne and especially when to switch to Babais enumeration help description.  */
  enum enum_dComp dComp_arg;	/**< @brief how to compute d Sequence for LP's enumeration (default='success').  */
  char * dComp_orig;	/**< @brief how to compute d Sequence for LP's enumeration original value given at command line.  */
  const char *dComp_help; /**< @brief how to compute d Sequence for LP's enumeration help description.  */
  enum enum_rComp rComp_arg;	/**< @brief how to compute R Sequence for Length Pruning (default='length').  */
  char * rComp_orig;	/**< @brief how to compute R Sequence for Length Pruning original value given at command line.  */
  const char *rComp_help; /**< @brief how to compute R Sequence for Length Pruning help description.  */
  double factor_arg;	/**< @brief controls the number of iterations done during enumeration (default='1.5').  */
  char * factor_orig;	/**< @brief controls the number of iterations done during enumeration original value given at command line.  */
  const char *factor_help; /**< @brief controls the number of iterations done during enumeration help description.  */
  double factor_bin_arg;	/**< @brief controls the number of iterations done when using binary secret d sequences (default='1.0').  */
  char * factor_bin_orig;	/**< @brief controls the number of iterations done when using binary secret d sequences original value given at command line.  */
  const char *factor_bin_help; /**< @brief controls the number of iterations done when using binary secret d sequences help description.  */
  long factor_lvl_arg;	/**< @brief used to balance short running threads by increasing the overall number of threads (default='8').  */
  char * factor_lvl_orig;	/**< @brief used to balance short running threads by increasing the overall number of threads original value given at command line.  */
  const char *factor_lvl_help; /**< @brief used to balance short running threads by increasing the overall number of threads help description.  */
  double delta_arg;	/**< @brief delta of BKZ reduction (default='0.99').  */
  char * delta_orig;	/**< @brief delta of BKZ reduction original value given at command line.  */
  const char *delta_help; /**< @brief delta of BKZ reduction help description.  */
  int parallel_flag;	/**< @brief run parallel implementations of enumeration algorithms (default=off).  */
  const char *parallel_help; /**< @brief run parallel implementations of enumeration algorithms help description.  */
  int n_threads_arg;	/**< @brief number of threads to use in parallel implementations (default='0').  */
  char * n_threads_orig;	/**< @brief number of threads to use in parallel implementations original value given at command line.  */
  const char *n_threads_help; /**< @brief number of threads to use in parallel implementations help description.  */
  char * ifile_arg;	/**< @brief input file prefix, \"_{matrix,vector}.dat\" will be appended. (default='samples').  */
  char * ifile_orig;	/**< @brief input file prefix, \"_{matrix,vector}.dat\" will be appended. original value given at command line.  */
  const char *ifile_help; /**< @brief input file prefix, \"_{matrix,vector}.dat\" will be appended. help description.  */
  char * ofile_arg;	/**< @brief output file (default='solution_vector.dat').  */
  char * ofile_orig;	/**< @brief output file original value given at command line.  */
  const char *ofile_help; /**< @brief output file help description.  */
  int binary_secret_flag;	/**< @brief generate secrect from binary uniform distribution {0, 1} (default=off).  */
  const char *binary_secret_help; /**< @brief generate secrect from binary uniform distribution {0, 1} help description.  */
  int binary_a_flag;	/**< @brief generate A from binary uniform distribution {0, 1} (default=off).  */
  const char *binary_a_help; /**< @brief generate A from binary uniform distribution {0, 1} help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int dimension_given ;	/**< @brief Whether dimension was given.  */
  unsigned int modulus_given ;	/**< @brief Whether modulus was given.  */
  unsigned int sigma_given ;	/**< @brief Whether sigma was given.  */
  unsigned int beta_given ;	/**< @brief Whether beta was given.  */
  unsigned int enumeration_given ;	/**< @brief Whether enumeration was given.  */
  unsigned int babaiBound_given ;	/**< @brief Whether babaiBound was given.  */
  unsigned int dComp_given ;	/**< @brief Whether dComp was given.  */
  unsigned int rComp_given ;	/**< @brief Whether rComp was given.  */
  unsigned int factor_given ;	/**< @brief Whether factor was given.  */
  unsigned int factor_bin_given ;	/**< @brief Whether factor_bin was given.  */
  unsigned int factor_lvl_given ;	/**< @brief Whether factor_lvl was given.  */
  unsigned int delta_given ;	/**< @brief Whether delta was given.  */
  unsigned int parallel_given ;	/**< @brief Whether parallel was given.  */
  unsigned int n_threads_given ;	/**< @brief Whether n-threads was given.  */
  unsigned int ifile_given ;	/**< @brief Whether ifile was given.  */
  unsigned int ofile_given ;	/**< @brief Whether ofile was given.  */
  unsigned int binary_secret_given ;	/**< @brief Whether binary_secret was given.  */
  unsigned int binary_a_given ;	/**< @brief Whether binary_a was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern const char *cmdline_parser_enumeration_values[];  /**< @brief Possible values for enumeration. */
extern const char *cmdline_parser_dComp_values[];  /**< @brief Possible values for dComp. */
extern const char *cmdline_parser_rComp_values[];  /**< @brief Possible values for rComp. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_ENUMERATION_H */
