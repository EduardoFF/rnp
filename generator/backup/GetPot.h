//  -*- c++ -*-  vim: set syntax=cpp:
//  GetPot Version $$Version$$                                    $$Date$$
//
//  WEBSITE: http://getpot.sourceforge.net
//
//  NOTE: The LPGL License for this library is only valid in case that 
//        it is not used for the production or development of applications
//        dedicated to military industry. This is what the author calls
//        the 'unofficial peace version of the LPGL'.
//
//  This library is  free software; you can redistribute  it and/or modify
//  it  under  the terms  of  the GNU  Lesser  General  Public License  as
//  published by the  Free Software Foundation; either version  2.1 of the
//  License, or (at your option) any later version.
//
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//
//  You  should have  received a  copy of  the GNU  Lesser  General Public
//  License along  with this library; if  not, write to  the Free Software
//  Foundation, Inc.,  59 Temple Place,  Suite 330, Boston,  MA 02111-1307
//  USA
//
//  (C) 2001-2007 Frank R. Schaefer <fschaef@users.sf.net>
//==========================================================================

#ifndef __include_guard_GETPOT_H__
#define __include_guard_GETPOT_H__

#if defined(WIN32) || defined(SOLARIS_RAW) || (__GNUC__ == 2) || defined(__HP_aCC)
#define strtok_r(a, b, c) strtok(a, b)
#endif // WINDOWS or SOLARIS or gcc 2.* or HP aCC

extern "C" {
//   leave the 'extern C' to make it 100% sure to work -
//   expecially with older distributions of header files.
#ifndef WIN32
// this is necessary (depending on OS)
#include <ctype.h>
#endif
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
}
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include <fstream>
#include <iostream> // not every compiler distribution includes <iostream> 
//                  // with <fstream>

typedef  std::vector<std::string>  STRING_VECTOR;

#define victorate(TYPE, VARIABLE, ITERATOR)                        \
  std::vector<TYPE>::const_iterator ITERATOR = (VARIABLE).begin(); \
  for(; (ITERATOR) != (VARIABLE).end(); (ITERATOR)++)


class GetPot {
    //--------
     void __basic_initialization();
public:
    // (*) constructors, destructor, assignment operator -----------------------
     GetPot();
     GetPot(const GetPot&);
     GetPot(const int argc_, char** argv_, 
                  const char* FieldSeparator=0x0);
     GetPot(const char* FileName, 
                  const char* CommentStart=0x0, const char* CommentEnd=0x0,
                  const char* FieldSeparator=0x0);
     ~GetPot();
     GetPot& operator=(const GetPot&);


    // (*) absorbing contents of another GetPot object
     void            absorb(const GetPot& That);
    //     -- for ufo detection: recording requested arguments, options etc.
     void            clear_requests();
    inline void            disable_request_recording() { __request_recording_f = false; }
    inline void            enable_request_recording()  { __request_recording_f = true; }

    // (*) direct access to command line arguments -----------------------------
     const std::string    operator[](unsigned Idx) const;
     int                  get(unsigned Idx, int           Default) const;
     double               get(unsigned Idx, const double& Default) const;
     const std::string    get(unsigned Idx, const char*   Default) const;
     unsigned             size() const;

    // (*) flags ---------------------------------------------------------------
     bool            options_contain(const char* FlagList) const;
     bool            argument_contains(unsigned Idx, const char* FlagList) const;

    // (*) variables -----------------------------------------------------------
    //     -- scalar values
     int                operator()(const char* VarName, int           Default) const;
     double             operator()(const char* VarName, const double& Default) const;
     const std::string  operator()(const char* VarName, const char*   Default) const;
    //     -- vectors
     int                operator()(const char* VarName, int Default, unsigned Idx) const;
     double             operator()(const char* VarName, const double& Default, unsigned Idx) const;
     const std::string  operator()(const char* VarName, const char* Default, unsigned Idx) const;

    //     -- setting variables
    //                  i) from outside of GetPot (considering prefix etc.)
    //                  ii) from inside, use '__set_variable()' below
     void            set(const char* VarName, const char* Value, const bool Requested = true);
     void            set(const char* VarName, const double& Value, const bool Requested = true);
     void            set(const char* VarName, const int Value, const bool Requested = true);
    
     unsigned        vector_variable_size(const char* VarName) const;
     STRING_VECTOR   get_variable_names() const;
     STRING_VECTOR   get_section_names() const;
    

    // (*) cursor oriented functions -------------------------------------------
    inline void            set_prefix(const char* Prefix) { prefix = std::string(Prefix); }
    inline bool            search_failed() const { return search_failed_f; }

    //     -- enable/disable search for an option in loop
    inline void            disable_loop() { search_loop_f = false; }
    inline void            enable_loop()  { search_loop_f = true; }

    //     -- reset cursor to position '1'
     void            reset_cursor();
     void            init_multiple_occurrence();

    //     -- search for a certain option and set cursor to position
    bool            search(const char* option);
    bool            search(unsigned No, const char* P, ...);
    //     -- get argument at cursor++
     int                next(int           Default);
     double             next(const double& Default);
     const std::string  next(const char*   Default);
    //     -- search for option and get argument at cursor++
     int                follow(int           Default, const char* Option);
     double             follow(const double& Default, const char* Option);
     const std::string  follow(const char*   Default, const char* Option);
    //     -- search for one of the given options and get argument that follows it
     int                follow(int           Default, unsigned No, const char* Option, ...);
     double             follow(const double& Default, unsigned No, const char* Option, ...);
     const std::string  follow(const char*   Default, unsigned No, const char* Option, ...);
    //     -- lists of nominuses following an option
     std::vector<std::string> nominus_followers(const char* Option);
     std::vector<std::string> nominus_followers(unsigned No, ...);

    //     -- directly followed arguments
     int                direct_follow(int           Default, const char* Option);
     double             direct_follow(const double& Default, const char* Option);
     const std::string  direct_follow(const char*   Default, const char* Option);

     std::vector<std::string>  string_tails(const char* StartString);
     std::vector<int>          int_tails(const char* StartString, const int Default = 1);
     std::vector<double>       double_tails(const char* StartString, const double Default = 1.0);

    // (*) nominus arguments ---------------------------------------------------
     STRING_VECTOR   nominus_vector() const;
     unsigned        nominus_size() const  { return static_cast<unsigned int>(idx_nominus.size()); }
     std::string     next_nominus();

    // (*) unidentified flying objects -----------------------------------------
     STRING_VECTOR   unidentified_arguments(unsigned Number, const char* Known, ...) const;
     STRING_VECTOR   unidentified_arguments(const STRING_VECTOR& Knowns) const;
     STRING_VECTOR   unidentified_arguments() const;

     STRING_VECTOR   unidentified_options(unsigned Number, const char* Known, ...) const;
     STRING_VECTOR   unidentified_options(const STRING_VECTOR& Knowns) const;
     STRING_VECTOR   unidentified_options() const;

     std::string     unidentified_flags(const char* Known,
                                             int ArgumentNumber /* =-1 */) const;

     STRING_VECTOR   unidentified_variables(unsigned Number, const char* Known, ...) const;
     STRING_VECTOR   unidentified_variables(const STRING_VECTOR& Knowns) const;
     STRING_VECTOR   unidentified_variables() const;

     STRING_VECTOR   unidentified_sections(unsigned Number, const char* Known, ...) const;
     STRING_VECTOR   unidentified_sections(const STRING_VECTOR& Knowns) const;
     STRING_VECTOR   unidentified_sections() const;

     STRING_VECTOR   unidentified_nominuses(unsigned Number, const char* Known, ...) const;
     STRING_VECTOR   unidentified_nominuses(const STRING_VECTOR& Knowns) const;
     STRING_VECTOR   unidentified_nominuses() const;

    // (*) output --------------------------------------------------------------
     int print() const;

private:
    // (*) Type Declaration ----------------------------------------------------
    struct variable {
        //-----------
        // Variable to be specified on the command line or in input files.
        // (i.e. of the form var='12 312 341')

        // -- constructors, destructors, assignment operator
        variable();
        variable(const variable&);
        variable(const char* Name, const char* Value, const char* FieldSeparator);
        ~variable();
        variable& operator=(const variable& That);

        void      take(const char* Value, const char* FieldSeparator);

        // -- get a specific element in the string vector
        //    (return 0 if not present)
        const std::string*  get_element(unsigned Idx) const;

        // -- data memebers
        std::string       name;      // identifier of variable
        STRING_VECTOR     value;     // value of variable stored in vector
        std::string       original;  // value of variable as given on command line
    };

    // (*) member variables --------------------------------------------------------------
    std::string           prefix;          // prefix automatically added in queries
    std::string           section;         // (for dollar bracket parsing)
    STRING_VECTOR         section_list;    // list of all parsed sections
    //     -- argument vector
    STRING_VECTOR         argv;            // vector of command line arguments stored as strings
    unsigned              cursor;          // cursor for argv
    bool                  search_loop_f;   // shall search start at beginning after
    //                                     // reaching end of arg array ?
    bool                  search_failed_f; // flag indicating a failed search() operation
    //                                     // (e.g. next() functions react with 'missed')

    //     --  nominus vector
    int                   nominus_cursor;  // cursor for nominus_pointers
    std::vector<unsigned> idx_nominus;     // indecies of 'no minus' arguments

    //     -- variables
    //       (arguments of the form "variable=value")
    std::vector<variable> variables;
    
    //     -- comment delimiters
    std::string           _comment_start;
    std::string           _comment_end;

    //     -- field separator (separating elements of a vector)
    std::string           _field_separator;

    //     -- some functions return a char pointer to a temporarily existing string
    //        this container makes them 'available' until the getpot object is destroyed.
    std::vector<char*>    __internal_string_container;

    //     -- keeping track about arguments that are requested, so that the UFO detection
    //        can be simplified
    STRING_VECTOR   _requested_arguments;
    STRING_VECTOR   _requested_variables;
    STRING_VECTOR   _requested_sections;

    bool            __request_recording_f;   // speed: request recording can be turned off

    //     -- if an argument is requested record it and the 'tag' the section branch to which 
    //        it belongs. Caution: both functions mark the sections as 'tagged'.
    void                      __record_argument_request(const std::string& Arg);
    void                      __record_variable_request(const std::string& Arg);

    // (*) helper functions ----------------------------------------------------
    //                  set variable from inside GetPot (no prefix considered)
     void               __set_variable(const char* VarName, const char* Value);

    //     -- produce three basic data vectors:
    //          - argument vector
    //          - nominus vector
    //          - variable dictionary
     void               __parse_argument_vector(const STRING_VECTOR& ARGV);

    //     -- helpers for argument list processing
    //        * search for a variable in 'variables' array
     const variable*    __find_variable(const char*) const;
    //        * support finding directly followed arguments
     const char*        __match_starting_string(const char* StartString);
    //        * support search for flags in a specific argument
     bool               __check_flags(const std::string& Str, const char* FlagList) const;
    //        * type conversion if possible
     int                __convert_to_type(const std::string& String, int Default) const;
     double             __convert_to_type(const std::string& String, double Default) const;
    //        * prefix extraction
    const std::string         __get_remaining_string(const std::string& String, 
                                                     const std::string& Start) const;
    //        * search for a specific string
     bool               __search_string_vector(const STRING_VECTOR& Vec,
                                                     const std::string& Str) const;

    //     -- helpers to parse input file
    //        create an argument vector based on data found in an input file, i.e.:
    //           1) delete comments (in between '_comment_start' '_comment_end')
    //           2) contract assignment expressions, such as
    //                   my-variable   =    '007 J. B.'
    //             into
    //                   my-variable='007 J. B.'
    //           3) interprete sections like '[../my-section]' etc.
     void               __skip_whitespace(std::istream& istr);
     const std::string  __get_next_token(std::istream& istr);
     const std::string  __get_string(std::istream& istr);
     const std::string  __get_until_closing_bracket(std::istream& istr);

     STRING_VECTOR      __read_in_stream(std::istream& istr);
     STRING_VECTOR      __read_in_file(const char* FileName);
     std::string        __process_section_label(const std::string& Section,
                                                      STRING_VECTOR& section_stack);

    //      -- dollar bracket expressions
    std::string               __DBE_expand_string(const std::string str);
    std::string               __DBE_expand(const std::string str);
    const GetPot::variable*   __DBE_get_variable(const std::string str);
    STRING_VECTOR             __DBE_get_expr_list(const std::string str, const unsigned ExpectedNumber);

    std::string  __double2string(const double& Value) const {
        // -- converts a double integer into a string
        char* tmp = new char[128];
#ifndef WIN32
        snprintf(tmp, (int)sizeof(char)*128, "%e", Value);
#else
        _snprintf(tmp, sizeof(char)*128, "%e", Value);
#endif
        std::string result(tmp);
        delete [] tmp;
        return result;
    }

    std::string  __int2string(const int& Value) const {
        // -- converts an integer into a string
        char* tmp = new char[128];
#ifndef WIN32
        snprintf(tmp, (int)sizeof(char)*128, "%i", Value);
#else
        _snprintf(tmp, sizeof(char)*128, "%i", Value);
#endif
        std::string result(tmp);
        delete [] tmp;
        return result;
    }

    STRING_VECTOR __get_section_tree(const std::string& FullPath) {
        // -- cuts a variable name into a tree of sub-sections. this is requested for recording
        //    requested sections when dealing with 'ufo' detection.
        STRING_VECTOR   result;
        const char* Start = FullPath.c_str();

        for(char *p = (char*)Start; *p ; p++) {
            if( *p == '/' ) { 
                *p = '\0';  // set terminating zero for convinience
                const std::string Section = Start;
                *p = '/';   // reset slash at place
                result.push_back(Section);
            }
        }

        return result;
    }
};

#endif // __include_guard_GETPOT_H__



