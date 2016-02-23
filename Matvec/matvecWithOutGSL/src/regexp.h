#ifndef MATVEC_REGEXP_H
#define MATVEC_REGEXP_H

#include <string>
#include <vector>
#include <stdexcept>
#include <pcre.h>

#include "exception.h";
#include "util.h"
namespace matvec {
/**
 * A c++ class wrapper for package pcre, a Perl compatible reggular expression library.
 * You must have pcre-3.4 or better installed in your system. Add -lpcrc to your compiling
 * options. Note that you don't need -lmatvec to access PCRE.
 *
 * PCRE was developped by Philip Hazel<ph10@cam.ac.uk> at University of Cambridge.
 * The following is what he wrote:
 *
 *    The PCRE library is a set of functions that implement  regu-
 *    lar  expression  pattern  matching using the same syntax and
 *    semantics as Perl  5,  with  just  a  few  differences.
 *    The  current  implementation  corresponds  to  Perl
 *    5.005, with some additional features  from  later  versions.
 *
 *     A regular expression is a pattern that is matched against  a
 *     subject string from left to right. Most characters stand for
 *     themselves in a pattern, and match the corresponding charac-
 *     ters in the subject. As a trivial example, the pattern
 *
 *       The quick brown fox
 *
 *     matches a portion of a subject string that is  identical  to
 *     itself.  The  power  of  regular  expressions comes from the
 *     ability to include alternatives and repetitions in the  pat-
 *     tern.  These  are encoded in the pattern by the use of meta-
 *     characters, which do not stand for  themselves  but  instead
 *     are interpreted in some special way.
 *
 *     There are two different sets of meta-characters: those  that
 *     are  recognized anywhere in the pattern except within square
 *     brackets, and those that are recognized in square  brackets.
 *     Outside square brackets, the meta-characters are as follows:
 *
 * <UL>
 *   <LI> \      general escape character with several uses
 *   <LI> ^      assert start of  subject  (or  line,  in  multiline mode)
 *   <LI> $      assert end of subject (or line, in multiline mode)
 *                  match any character except newline (by default)
 *   <LI> [      start character class definition
 *   <LI> |      start of alternative branch
 *   <LI> (      start subpattern
 *   <LI> )      end subpattern
 *   <LI> ?      extends the meaning of (, also 0 or 1 quantifier, also quantifier minimizer
 *   <LI> *      0 or more quantifier
 *   <LI> +      1 or more quantifier
 *   <LI> {      start min/max quantifier
 * </UL>
 *
 *     Part of a pattern that is in square  brackets  is  called  a
 *     "character  class".  In  a  character  class  the only meta-
 *     characters are:
 * <UL>
 *   <LI>  \      general escape character
 *   <LI>  ^      negate the class, but only if the first character
 *   <LI>  -      indicates character range
 *   <LI>  ]      terminates the characte
 * </UL>
 *
 * Here are some examples:
 * <OL>
 *   <LI> "[12[:^digit:]]"  matches "1", "2", or any non-digit.
 *   <LI> "gilbert|sullivan" atches either "gilbert" or "sullivan".
 *   <LI> "(a(?i)b)c" matches  abc  and  aBc  and  no  other   strings
 *   <LI> "z{2,4}" matches "zz", "zzz", or "zzzz"
 *   <LI> "\d{8}" matches exactly 8 digits
 *   <LI> "*"    is equivalent to {0,}
 *   <LI> "+"    is equivalent to {1,}
 *   <LI> "?"    is equivalent to {0,1}
 *   <LI> "((?i)rah)\s+\1"  matches "rah rah" and "RAH RAH", but  not  "RAH  rah"
 *   <LI> "(a|b\1)+" matches any number of "a"s and also "aba", "ababbaa" etc.
 * </OL>
 * </CODE>
 */
class RegExp
{
public:
   enum // bits of opciones
   {
      anchored  = PCRE_ANCHORED,
      caseless  = PCRE_CASELESS,
      dollarend = PCRE_DOLLAR_ENDONLY,
      dotall    = PCRE_DOTALL,
      extended  = PCRE_EXTENDED,
      multiline = PCRE_MULTILINE,
      ungreedy  = PCRE_UNGREEDY
   };

   RegExp(unsigned opts = 0)                       { re = 0; _opts = opts;}             ///< Constructs a regexp
   RegExp(const std::string &s, unsigned opts = 0) { re = 0; _opts = opts; compile(s);} ///< Constructs a regexp
   RegExp(const char *s, unsigned opts = 0)        { re = 0; _opts = opts; compile(s);} ///< Constructs a regexp
   RegExp(const RegExp &r)                         { re = 0; _opts = 0; copy(r); }      ///< A copu constructor
   ~RegExp() { if ( re ) delete[] re; re = 0; }                                         ///< Destruts the object

   const RegExp & operator = (const std::string &s) { compile(s); return *this; }       ///< Assignment.
   const RegExp & operator = (const char *s)        { compile(s); return *this; }       ///< Assignment.
   const RegExp & operator = (const RegExp &r)      { copy(r); return *this; }          ///< Assignment.


   unsigned		             options() { return _opts; }                       ///< Returns the options
   void			             options(unsigned opts) { _opts = opts; }          ///< Set the options
   int                               find(const std::string &s, std::vector<std::string> *subs = 0, const unsigned offset = 0) const;
   int                               ngmatch(const std::string &s) const;
   std::string                       replace(const std::string &s,const std::string &rep,const unsigned offset = 0) const;
   std::vector<std::pair<int, int> > match(const std::string &s, unsigned offset = 0) const;
   std::vector<std::pair<int, int> > gmatch(const std::string &s) const;
   static std::string                substr(const std::string &s,const std::vector<std::pair<int, int> > &marks, unsigned index);
   static std::vector<std::string>   substr(const std::string &s,const std::vector<std::pair<int, int> > &marks);
   std::vector<std::string>          split(const std::string &s, bool emptyfields = true,unsigned maxfields = 0) const;
private:
   pcre *re;
   unsigned _opts;

   void compile(const std::string &s);
   void copy(const RegExp &r);
};

inline void RegExp::compile(const std::string &s)
{
   const char * errorptr;
   int erroroffset;

   if ( re ) delete[] re;
   re = pcre_compile(s.c_str(), _opts, &errorptr, &erroroffset, 0);
   if ( !re ) throw exception(string(errorptr) + " at: " + s.substr(erroroffset));
}

inline void RegExp::copy(const RegExp &r)
{
   if (this == &r) return;
   size_t size;
   pcre_fullinfo(r.re, 0, PCRE_INFO_SIZE, &size);
   if ( re ) delete[] re;
   if(size>0){
     re = (pcre *) new char[size];
   }
   else {
     re = 0;
   }
   check_ptr(re);
   memcpy(re, r.re, size);
   _opts = r._opts;
   return;
}

/**
 * Matches it against the compiled regular expression.
 *
 * It takes a string and an optional starting offset, whose default value is 0.
 * An exception is thrown if the Regexp is uninitialized. It returns a vector
 * of pair<int,int>. If the returned vector is empty, the string did not match.
 *  If the returned vector (let's call it v) is not empty,
 *  then the v[0] pair contains the offsets
 * of the first and the last-plus-one characters in the match. The rest of
 * the elements, v[i], contain the same information for the captured substrings.
 * If a certain subpattern in the expression did not participate in the match,
 * the corresponding vector element will contain the pair (-1, -1).
 */
inline std::vector<std::pair<int,int> > RegExp::match(const std::string &s, unsigned offset) const
{
   if ( !re ) throw exception("match on uninitialized expression");

   size_t msize;
   pcre_fullinfo(re, 0, PCRE_INFO_CAPTURECOUNT, &msize);
   msize = 3*(msize+1);
   int *m = new int[msize];
   check_ptr(m);
   int result;
   vector<pair<int,int> > marks;

   result = pcre_exec(re, 0, s.c_str(), s.length(), offset, 0, m, msize);
   for ( int i = 0, *p = m ; i < result ; i++, p+=2 ) marks.push_back(make_pair(p[0], p[1]));
   delete[] m;
   return marks;
}

/**
 * Finds all matchings of the regular expression
 * in the string and returns a vector of pair<int,int> with one element
 * for each match. Substrings for each match are not reported.
 * The strings corresponding to the individual matches
 * can be retrieved using the substr() functions.
 */
inline std::vector<std::pair<int,int> > RegExp::gmatch(const std::string &s) const
{
   if ( !re ) throw exception("gmatch on uninitialized expression");

   int m[3];
   vector<pair<int,int> > marks;

   const char * str = s.c_str();
   unsigned offset = 0, len = s.length();
   while ( offset < len && pcre_exec(re, 0, str, len, offset, 0, m, 3) >= 0 ) {
      marks.push_back(make_pair(m[0], m[1]));
      offset = m[1];
   }
   return marks;
}

/**
 * It returns the corresponding substring.
 * The returned substring is corresponding to
 *  the matched string, the vector of pair<int,int>
 * returned by match() and an index.
 *
 */
inline std::string RegExp::substr(const std::string &s, const std::vector<std::pair<int,int> > &marks,unsigned index)
{
   if ( index >= marks.size() ) throw exception("bad substring index");

   int begin = marks[index].first;
   if ( begin == -1 ) return "";
   int end = marks[index].second;
   return s.substr(begin, end-begin);
}

/**
 * Returns a vector of strings containing all the substrings in the match.
 */
inline std::vector<std::string> RegExp::substr(const std::string & s, const std::vector<std::pair<int,int> > & marks)
{
   vector<string> v;
   unsigned size = marks.size();

   for ( unsigned i = 0 ; i < size ; i++ ) v.push_back(substr(s, marks, i));
   return v;
}

/**
 * Returns a vector of substrings splitted from string s.
 */
inline std::vector<std::string> RegExp::split(const std::string &s, bool emptyfields, unsigned maxfields) const
{
   if ( !re ) throw exception("split on uninitialized expression");
   vector<pair<int,int> > m = gmatch(s);
   vector<pair<int,int> > marks;

   int begin = 0, end;
   for ( int i = 0, nsep = m.size() ; i < nsep ; i++ ) {
      end = m[i].first;
      if ( emptyfields || end > begin )
         marks.push_back(make_pair(begin, end));
      begin = m[i].second;
   }
   end = s.length();
   if ( emptyfields || end > begin ) marks.push_back(make_pair(begin, end));
   unsigned nfields = marks.size();
   if ( maxfields && nfields > maxfields ) {
      marks[maxfields-1].second = marks[nfields-1].second;
      marks.erase(&marks[maxfields], marks.end());
   }
   return substr(s, marks);
}

/**
 * Try to find the regexp by returning the index of first occurrence.
 * It brings out a vector of all findings.

 */
inline int RegExp::find(const std::string &s, std::vector<std::string> *subs, const unsigned offset) const
{
   if ( !re ) throw exception("find on uninitialized expression");

   size_t msize;
   pcre_fullinfo(re, 0, PCRE_INFO_CAPTURECOUNT, &msize);
   msize = 3*(msize+1);
   int *m = new int[msize];
   check_ptr(m);
   int result = pcre_exec(re, 0, s.c_str(), s.length(), offset, 0, m, msize);
   int ret;
   if (result >= 0) {
      ret = m[0];
   } else {
      ret = result;
   }

   if (subs) {
      if (result > 0) {
         subs->reserve(result);
         for (int j,i = 0; i < result ; i++) {
             j = i*2;
             if (m[j] == -1) {
                subs->push_back("");
             } else {
                subs->push_back(s.substr(m[j], m[j+1] - m[j]));
             }
         }
      } else {
         subs->clear();
      }
   }
   delete[] m;
   return ret;
}

/**
 *  Returns the number of matches
 */
inline int RegExp::ngmatch(const std::string &s) const
{
   if ( !re ) throw exception("nmatch on uninitialized expression");

   int ret = 0;
   int m[3];
   const char *str = s.c_str();
   unsigned offset = 0, len = s.length();
   while ( offset < len && pcre_exec(re, 0, str, len, offset, 0, m, 3) >= 0 ) {
      ret++;
      offset = m[1];
   }
   return ret++;
}

/**
 * Replace occurrence with rep starting from offset.
 */
inline std::string RegExp::replace(const std::string &str,const std::string &rep,const unsigned offset) const
{
   if ( !re ) throw exception("replace on uninitialized expression");

   size_t msize;
   pcre_fullinfo(re, 0, PCRE_INFO_CAPTURECOUNT, &msize);
   msize = 3*(msize+1);
   int *m = new int[msize];
   check_ptr(m);
   int i,j,k,len,nmat;
   len = str.length();
   if ((nmat = pcre_exec(re, 0, str.c_str(), len, offset, 0, m, msize)) < 0) return str;

   string ret = rep;
   string::size_type begidx;
   begidx = ret.find_first_of("$",0);
   while (begidx != string::npos) {
      k = 1;
      if (isdigit(ret.at(begidx+1))) {
         i = ret[begidx + 1]- 48;
         if (i && i < nmat) {
            j = i*2;
            k = m[j+1] - m[j];
            ret.replace(begidx,2,str.substr(m[j],k));
         }
      }
      begidx = ret.find_first_of("$",begidx + k);
   }

   if (m[1] < len) ret.append(str.substr(m[1],len - m[1]));
   if (m[0] > 0) ret = str.substr(0,m[0])  + ret;
   return ret;
}
} //////// end of namespace matvec
#endif
