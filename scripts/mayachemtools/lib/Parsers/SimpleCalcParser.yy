%{
package Parsers::SimpleCalcParser;
#
# File: SimpleCalcParser.yy
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2017 Manish Sud. All rights reserved.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#
#
# A WORD TO THE WISE:
#
# The parser package and token table files, SimpleCalcParser.pm and
#  SimpleCalcParser.tab.ph, are automatically generated from parser grammar
# definition file, SimpleCalcParser.yy, using byacc available through perl-byacc1.8
# modified with perl5-byacc-patches-0.5 for generation of object oriented parser:
#
#    byacc -l -P -d -b SimpleCalcParser SimpleCalcParser.yy
#    mv SimpleCalcParser.tab.pl SimpleCalcParser.pm
#

use Carp;

# Setup a hash map for mapping of words/letters to values...
%LetterToValueMap = ();

%}

%start list

%token NUMBER LETTER

%left '+' '-'
%left '*' '/'
%left '%'

%%

list    :
        |       list  stat  '\n'
                        {  $$ = $2; }
        |       list  error  '\n'
                        { $p->yyerrok; $p->yyclearin; }
        ;


stat    :       expr
                        {  $ExprOut = sprintf "%5i", $1; print "$ExprOut\n"; $$ = $1; }
        |       LETTER  '='  expr
                        { $LetterToValueMap{$1} = $3; }
        ;

expr    :       '('  expr  ')'
                        { $$ = $2; }
        |       expr '+' expr
                        { $$ = $1 + $3; }
        |       expr '-' expr
                        { $$ = $1 - $3; }
        |       expr '*' expr
                        { $$ = $1 * $3; }
        |       expr '/' expr
                        { $$ = $1 / $3; }
        |       expr '%' expr
                        { $$ = $1 % $3; }
        |       NUMBER
        |       LETTER
                        {
                          if (exists $LetterToValueMap{$1}) {
                            $$ = $LetterToValueMap{$1};
                          }
                          else {
                            $Letter = $1;
                            print "Undefined variable $Letter encountered by SimpleCalcParser; Value set to 0\n";
                            $$ = 0;
                          }
                        }
        ;

%%

# yyerror function supplied to parser along with a lexer during initialization of
# the parser...
#
sub yyerror {
    my ($msg) = @_;
    print "yyerror: $msg...\n";
}

1;

__END__

=head1 NAME

Parsers::SimpleCalcParser

=head1 SYNOPSIS

use Parsers::SimpleCalcParser ;

use Parsers::SimpleCalcParser qw(:all);

=head1 DESCRIPTION

B<Parsers::SimpleCalcParser> class provides the following methods:

new, yyclearin, yyerrok, yyerror, yyparse

B<Parsers::SimpleCalcParse.yy> parser grammer definition file implements a simple
calculator and is provided to highlight usage of lexer capability available through
B<Parsers::SimpleCalcYYLexer>, which in turn uses B<Parsers::YYLexer> and
B<Parsers::Lexer> classes to procide underlying lexer functionality.

The parser package and token table files, B<Parsers::SimpleCalcParser.pm> and
B<SimpleCalcParser.tab.ph>, are automatically generated from parser grammar definition
file, B<Parsers::SimpleCalcParser.yy>, using byacc available through perl-byacc1.8 modified
with perl5-byacc-patches-0.5 for generation of object oriented parser:

    byacc -l -P -d -b SimpleCalcParser SimpleCalcParser.yy
    mv SimpleCalcParser.tab.pl SimpleCalcParser.pm

=head2 METHODS

=over 4

=item B<new>

    $SimpleCalcParser = new Parsers::SimpleCalcParser($YYLex,
                                \&Parsers::SimpleCalcParser::yyerror);
    $SimpleCalcParser = new Parsers::SimpleCalcParser($YYLex,
                                \&Parsers::SimpleCalcParser::yyerror, $Debug);

Using specified I<YYLex> I<YYError> functions, B<new> method generates a new
B<SimpleCalcParser> and returns a reference to newly created B<SimpleCalcYYParser> object.

Examples:

    # Input string...
    $InputText = "3 + 4 +6\nx=3\ny=5\nx+y\nx+z\n";
    $YYLexer = new Parsers::SimpleCalcYYLexer($InputText);
    $YYLex = $YYLexer->GetYYLex();

    $Debug = 0;
    $SimpleCalcParser = new Parsers::SimpleCalcParser($YYLex,
                               \&Parsers::SimpleCalcParser::yyerror, $Debug);
    $Value = $SimpleCalcParser->yyparse();
    print "Value = " . (defined($Value) ? "$Value" : "Undefined") . "\n";

    # Input file...
    $InputFile = "TestSimpleCalcParser.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";

    $YYLexer = new Parsers::SimpleCalcYYLexer(\*INPUTFILE);
    $YYLex = $YYLexer->GetYYLex();

    $Debug = 0;
    $SimpleCalcParser = new Parsers::SimpleCalcParser($YYLex,
                               \&Parsers::SimpleCalcParser::yyerror, $Debug);
    $Value = $SimpleCalcParser->yyparse();
    print "Value = " . (defined($Value) ? "$Value" : "Undefined") . "\n";

    close INPUTFILE;

    # Input iterator...
    $InputFile = "TestSimpleCalcParser.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $InputIterator = sub { return <INPUTFILE>; };

    $YYLexer = new Parsers::SimpleCalcYYLexer($InputIterator);
    $YYLex = $YYLexer->GetYYLex();

    $Debug = 0;
    $SimpleCalcParser = new Parsers::SimpleCalcParser($YYLex,
                               \&Parsers::SimpleCalcParser::yyerror, $Debug);
    $Value = $SimpleCalcParser->yyparse();
    print "Value = " . (defined($Value) ? "$Value" : "Undefined") . "\n";

    close INPUTFILE;

=item B<yyclearin>

    $SimpleCalcParser->yyclearin();

B<yyclearin> method clears any previous look-ahead token after encountering a syntax error
during parsing. It can be used after B<yyerrok> in a grammer rule with the reserved word
B<error>.

=item B<yyerrok>

    $SimpleCalcParser->yyerrok();

B<yyerrok> method is used with the reserved word B<error> in grammer rule to indcate
error recovery is complete after encountering a syntax error during parsing.

=item B<yyerror>

    $SimpleCalcParser->yyerror();

B<yyerror> function is provided for the caller to use during initialization of a parser. It
is used by B<yyparse> to print any error messages encountered during parsing of the
input.

=item B<yyparse>

    $Value = $SimpleCalcParser->yyparse();

Returns I<Value> after parsing all the input from a input stream using specified
grammer rules.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Lexer.pm, YYLexer.pm, SimpleCalcYYLexer.pm

=head1 COPYRIGHT

Copyright (C) 2017 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut

