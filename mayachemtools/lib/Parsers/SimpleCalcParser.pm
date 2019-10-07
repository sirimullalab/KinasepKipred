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

$NUMBER=257;
$LETTER=258;
$YYERRCODE=256;
@yylhs = (                                               -1,
    0,    0,    0,    1,    1,    2,    2,    2,    2,    2,
    2,    2,    2,
);
@yylen = (                                                2,
    0,    3,    3,    1,    3,    3,    3,    3,    3,    3,
    3,    1,    1,
);
@yydefred = (                                             1,
    0,    0,   12,    0,    0,    0,    0,    3,    0,   13,
    0,    2,    0,    0,    0,    0,    0,    0,    6,    0,
    0,    0,    0,   11,
);
@yydgoto = (                                              1,
    6,    7,
);
@yysindex = (                                             0,
  -40,   -7,    0,  -57,  -38,   -5,  -18,    0,  -38,    0,
  -31,    0,  -38,  -38,  -38,  -38,  -38,  -18,    0,  -16,
  -16,  -30,  -30,    0,
);
@yyrindex = (                                             0,
    0,    0,    0,   -9,    0,    0,   -1,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    3,    0,    8,
   13,   -2,    5,    0,
);
@yygindex = (                                             0,
    0,   50,
);
$YYTABLESIZE=220;
@yytable = (                                              5,
   13,    5,    8,    9,   12,   17,   17,    9,    4,   19,
   15,   13,    5,   14,   10,   16,    0,    7,   17,    0,
   17,    0,    8,   15,   13,   15,   14,   13,   16,    0,
   16,    0,   13,   13,    0,   13,    0,   13,    9,    9,
    9,    0,    9,    0,    9,   10,   10,   10,    7,   10,
    7,   10,    7,    8,   11,    8,    0,    8,   18,    0,
    0,    0,   20,   21,   22,   23,   24,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    2,    3,    4,    3,   10,
);
@yycheck = (                                             40,
   10,   40,   10,   61,   10,   37,   37,   10,   10,   41,
   42,   43,   10,   45,   10,   47,   -1,   10,   37,   -1,
   37,   -1,   10,   42,   43,   42,   45,   37,   47,   -1,
   47,   -1,   42,   43,   -1,   45,   -1,   47,   41,   42,
   43,   -1,   45,   -1,   47,   41,   42,   43,   41,   45,
   43,   47,   45,   41,    5,   43,   -1,   45,    9,   -1,
   -1,   -1,   13,   14,   15,   16,   17,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,  256,  257,  258,  257,  258,
);
$YYFINAL=1;
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
$YYMAXTOKEN=258;
#if YYDEBUG
@yyname = (
"end-of-file",'','','','','','','','','',"'\\n'",'','','','','','','','','','','','','','','','','','','','',
'','','','','','',"'%'",'','',"'('","')'","'*'","'+'",'',"'-'",'',"'/'",'','','','','','','','','',
'','','','',"'='",'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','',"NUMBER","LETTER",
);
@yyrule = (
"\$accept : list",
"list :",
"list : list stat '\\n'",
"list : list error '\\n'",
"stat : expr",
"stat : LETTER '=' expr",
"expr : '(' expr ')'",
"expr : expr '+' expr",
"expr : expr '-' expr",
"expr : expr '*' expr",
"expr : expr '/' expr",
"expr : expr '%' expr",
"expr : NUMBER",
"expr : LETTER",
);
#endif

sub yyclearin { $_[0]->{'yychar'} = -1; }

sub yyerrok { $_[0]->{'yyerrflag'} = 0; }

sub new {
  my $p = {'yylex' => $_[1], 'yyerror' => $_[2], 'yydebug' => $_[3]};
  bless $p, $_[0];
}

sub YYERROR { ++$_[0]->{'yynerrs'}; $_[0]->yy_err_recover; }

sub yy_err_recover {
  #
  # msud@san.rr.com:
  #
  # Turn off "exiting" warning to suppress the following warning at "next yyloop":
  #
  # Exiting subroutine via next at <LineNum>
  #
  # The code does work as expected with or without turning off the warning.
  # This method is invoked in yyparse method directly or indirectly in another
  # method and Perl compilers ends up finding "yyloop" as the nearst enclosure
  # label.
  #
  no warnings qw(exiting);

  my ($p) = @_;
  if ($p->{'yyerrflag'} < 3)
  {
    $p->{'yyerrflag'} = 3;
    while (1)
    {
      if (($p->{'yyn'} = $yysindex[$p->{'yyss'}->[$p->{'yyssp'}]]) && 
          ($p->{'yyn'} += $YYERRCODE) >= 0 && 
          $p->{'yyn'} < @yycheck &&
          $yycheck[$p->{'yyn'}] == $YYERRCODE)
      {
        warn("yydebug: state " . 
                     $p->{'yyss'}->[$p->{'yyssp'}] . 
                     ", error recovery shifting to state" . 
                     $yytable[$p->{'yyn'}] . "\n") 
                       if $p->{'yydebug'};
        $p->{'yyss'}->[++$p->{'yyssp'}] = 
          $p->{'yystate'} = $yytable[$p->{'yyn'}];
        $p->{'yyvs'}->[++$p->{'yyvsp'}] = $p->{'yylval'};
        next yyloop;
      }
      else
      {
        warn("yydebug: error recovery discarding state ".
              $p->{'yyss'}->[$p->{'yyssp'}]. "\n") 
                if $p->{'yydebug'};
        return(undef) if $p->{'yyssp'} <= 0;
        --$p->{'yyssp'};
        --$p->{'yyvsp'};
      }
    }
  }
  else
  {
    return (undef) if $p->{'yychar'} == 0;
    if ($p->{'yydebug'})
    {
      $p->{'yys'} = '';
      if ($p->{'yychar'} <= $YYMAXTOKEN) { $p->{'yys'} = 
        $yyname[$p->{'yychar'}]; }
      if (!$p->{'yys'}) { $p->{'yys'} = 'illegal-symbol'; }
      warn("yydebug: state " . $p->{'yystate'} . 
                   ", error recovery discards " . 
                   "token " . $p->{'yychar'} . "(" . 
                   $p->{'yys'} . ")\n");
    }
    $p->{'yychar'} = -1;
    next yyloop;
  }
0;
} # yy_err_recover

sub yyparse {
  my ($p, $s) = @_;
  if ($p->{'yys'} = $ENV{'YYDEBUG'})
  {
    $p->{'yydebug'} = int($1) if $p->{'yys'} =~ /^(\d)/;
  }

  $p->{'yynerrs'} = 0;
  $p->{'yyerrflag'} = 0;
  $p->{'yychar'} = (-1);

  $p->{'yyssp'} = 0;
  $p->{'yyvsp'} = 0;
  $p->{'yyss'}->[$p->{'yyssp'}] = $p->{'yystate'} = 0;

  yyloop: while(1)
  {
    yyreduce: {
      last yyreduce if ($p->{'yyn'} = $yydefred[$p->{'yystate'}]);
      if ($p->{'yychar'} < 0)
      {
        if ((($p->{'yychar'}, $p->{'yylval'}) = 
            &{$p->{'yylex'}}($s)) < 0) { $p->{'yychar'} = 0; }
        if ($p->{'yydebug'})
        {
          $p->{'yys'} = '';
          if ($p->{'yychar'} <= $#yyname) 
             { $p->{'yys'} = $yyname[$p->{'yychar'}]; }
          if (!$p->{'yys'}) { $p->{'yys'} = 'illegal-symbol'; };
          warn("yydebug: state " . $p->{'yystate'} . 
                       ", reading " . $p->{'yychar'} . " (" . 
                       $p->{'yys'} . ")\n");
        }
      }
      if (($p->{'yyn'} = $yysindex[$p->{'yystate'}]) && 
          ($p->{'yyn'} += $p->{'yychar'}) >= 0 && 
          $yycheck[$p->{'yyn'}] == $p->{'yychar'})
      {
        warn("yydebug: state " . $p->{'yystate'} . 
                     ", shifting to state " .
              $yytable[$p->{'yyn'}] . "\n") if $p->{'yydebug'};
        $p->{'yyss'}->[++$p->{'yyssp'}] = $p->{'yystate'} = 
          $yytable[$p->{'yyn'}];
        $p->{'yyvs'}->[++$p->{'yyvsp'}] = $p->{'yylval'};
        $p->{'yychar'} = (-1);
        --$p->{'yyerrflag'} if $p->{'yyerrflag'} > 0;
        next yyloop;
      }
      if (($p->{'yyn'} = $yyrindex[$p->{'yystate'}]) && 
          ($p->{'yyn'} += $p->{'yychar'}) >= 0 &&
          $yycheck[$p->{'yyn'}] == $p->{'yychar'})
      {
        $p->{'yyn'} = $yytable[$p->{'yyn'}];
        last yyreduce;
      }
      if (! $p->{'yyerrflag'}) {
        if ( (defined($EOI) && $p->{'yychar'} == $EOI) || ($p->{'yychar'} == 0) ) {
          &{$p->{'yyerror'}}("syntax error at or near the end of input text", $s);
        }
        else {
          &{$p->{'yyerror'}}("syntax error at or near input text: '$p->{'yylval'}'", $s);
        }
        ++$p->{'yynerrs'};
      }
      return(undef) if $p->yy_err_recover;
    } # yyreduce
    warn("yydebug: state " . $p->{'yystate'} . 
                 ", reducing by rule " . 
                 $p->{'yyn'} . " (" . $yyrule[$p->{'yyn'}] . 
                 ")\n") if $p->{'yydebug'};
    $p->{'yym'} = $yylen[$p->{'yyn'}];
    $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}+1-$p->{'yym'}];

    if ($p->{'yyn'} == 2) {
    {  $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-1]; }
    }

    if ($p->{'yyn'} == 3) {
    { $p->yyerrok; $p->yyclearin; }
    }

    if ($p->{'yyn'} == 4) {
    {  $ExprOut = sprintf "%5i", $p->{'yyvs'}->[$p->{'yyvsp'}-0]; print "$ExprOut\n"; $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 5) {
    { $LetterToValueMap{$p->{'yyvs'}->[$p->{'yyvsp'}-2]} = $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 6) {
    { $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-1]; }
    }

    if ($p->{'yyn'} == 7) {
    { $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-2] + $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 8) {
    { $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-2] - $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 9) {
    { $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-2] * $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 10) {
    { $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-2] / $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 11) {
    { $p->{'yyval'} = $p->{'yyvs'}->[$p->{'yyvsp'}-2] % $p->{'yyvs'}->[$p->{'yyvsp'}-0]; }
    }

    if ($p->{'yyn'} == 13) {
    {
                          if (exists $LetterToValueMap{$p->{'yyvs'}->[$p->{'yyvsp'}-0]}) {
                            $p->{'yyval'} = $LetterToValueMap{$p->{'yyvs'}->[$p->{'yyvsp'}-0]};
                          }
                          else {
                            $Letter = $p->{'yyvs'}->[$p->{'yyvsp'}-0];
                            print "Undefined variable $Letter encountered by SimpleCalcParser; Value set to 0\n";
                            $p->{'yyval'} = 0;
                          }
                        }
    }
    $p->{'yyssp'} -= $p->{'yym'};
    $p->{'yystate'} = $p->{'yyss'}->[$p->{'yyssp'}];
    $p->{'yyvsp'} -= $p->{'yym'};
    $p->{'yym'} = $yylhs[$p->{'yyn'}];
    if ($p->{'yystate'} == 0 && $p->{'yym'} == 0)
    {
      warn("yydebug: after reduction, shifting from state 0 ",
            "to state $YYFINAL\n") if $p->{'yydebug'};
      $p->{'yystate'} = $YYFINAL;
      $p->{'yyss'}->[++$p->{'yyssp'}] = $YYFINAL;
      $p->{'yyvs'}->[++$p->{'yyvsp'}] = $p->{'yyval'};
      if ($p->{'yychar'} < 0)
      {
        if ((($p->{'yychar'}, $p->{'yylval'}) = 
            &{$p->{'yylex'}}($s)) < 0) { $p->{'yychar'} = 0; }
        if ($p->{'yydebug'})
        {
          $p->{'yys'} = '';
          if ($p->{'yychar'} <= $#yyname) 
            { $p->{'yys'} = $yyname[$p->{'yychar'}]; }
          if (!$p->{'yys'}) { $p->{'yys'} = 'illegal-symbol'; }
          warn("yydebug: state $YYFINAL, reading " . 
               $p->{'yychar'} . " (" . $p->{'yys'} . ")\n");
        }
      }
      return ($p->{'yyvs'}->[1]) if $p->{'yychar'} == 0;
      next yyloop;
    }
    if (($p->{'yyn'} = $yygindex[$p->{'yym'}]) && 
        ($p->{'yyn'} += $p->{'yystate'}) >= 0 && 
        $p->{'yyn'} <= $#yycheck && 
        $yycheck[$p->{'yyn'}] == $p->{'yystate'})
    {
        $p->{'yystate'} = $yytable[$p->{'yyn'}];
    } else {
        $p->{'yystate'} = $yydgoto[$p->{'yym'}];
    }
    warn("yydebug: after reduction, shifting from state " . 
        $p->{'yyss'}->[$p->{'yyssp'}] . " to state " . 
        $p->{'yystate'} . "\n") if $p->{'yydebug'};
    $p->{'yyss'}[++$p->{'yyssp'}] = $p->{'yystate'};
    $p->{'yyvs'}[++$p->{'yyvsp'}] = $p->{'yyval'};
  } # yyloop
} # yyparse


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

