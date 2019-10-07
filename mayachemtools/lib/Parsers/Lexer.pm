package Parsers::Lexer;
#
# File: Lexer.pm
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
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

use strict;
use Carp;
use Exporter;
use Scalar::Util ();

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyLexer';

# Class constructor...
sub new {
  my($Class, $Input, @TokensSpec) = @_;

  # Initialize object...
  my $This = {};
  bless $This, ref($Class) || $Class;
  $This->_InitializeLexer();

  $This->_ValidateParametersAndGenerateLexer($Input, @TokensSpec);

  return $This;
}


# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object data...
#
sub _InitializeLexer {
  my($This) = @_;

  # Input parameter used by lexer to retrieve text to be lexed. Supported parameter types:
  #   . Reference to input iterator function
  #   . Reference to an open file handle
  #   . Text string
  #
  $This->{Input} = undef;

  # Type of input paramater determined using Perl ref function:
  #   . InputIterator - ref returns CODE
  #   . FileStream - ref return GLOB and fileno is valid
  #   . String - ref return an empty string
  #
  $This->{InputType} = '';

  # Tokens specifications supplied by the caller. It's an array containing references
  # to arrays with each containing TokenLabel and TokenMatchRegex pair along with
  # an option reference to code to be executed after a matched.
  #
  # For example:
  #
  # @LexerTokensSpec = (
  #			[ 'LETTER', qr/[a-zA-Z]/ ],
  #			[ 'NUMBER', qr/\d+/ ],
  #			[ 'SPACE', qr/[ ]*/, sub { my($This, $TokenLabel, $MatchedText) = @_; return ''; } ],
  #			[ 'NEWLINE', qr/(?:\r\n|\r|\n)/, sub { my($This, $TokenLabel, $MatchedText) = @_;  return "\n"; } ],
  #			[ 'CHAR', qr/[\.]/ ],
  #		       );
  #
  @{$This->{TokensSpec}} = ();

  # Refernce to chained lexer...
  $This->{ChainedLexer} = undef;

  return $This;
}

# Validate input parameters and generate a chained lexer...
#
sub _ValidateParametersAndGenerateLexer {
  my($This, $Input, @TokensSpec) = @_;

  #
  # Validate input to be lexed...
  if (!defined $Input) {
    croak "Error: ${ClassName}->new: Object can't be instantiated: Input is not defined. Supported values: a reference to input iterator function, a reference to an open file handle or a text string...";
  }
  $This->{Input} = $Input;

  # Check input parameter type...
  my($InputType);

  $InputType = ref $Input;
  if ($InputType =~ /CODE/i) {
    # Input iterator...
    $This->{InputType} = "InputIterator";
  }
  elsif ($InputType =~ /GLOB/i && defined fileno $Input) {
    # Input stream...
    $This->{InputType} = "FileStream";
  }
  elsif ($InputType) {
    # Perl ref function returns nonempty string for all other references...
    croak "Error: ${ClassName}->new: Object can't be instantiated: Invalid input parameter type specified. Supported parameter types: a reference to input iterator function, a reference to an open file handle or a text string...";
  }
  else {
    # Input string...
    $This->{InputType} = "String";
  }

  # Check tokens specifications...
  if (!@TokensSpec) {
    croak "Error: ${ClassName}->new: TokensSpec is not defined or the array doesn't contain any values. Supported values: a reference to an array containg token label, regular expression to match and an option reference to function to modify matched values...";
  }
  @{$This->{TokensSpec}} = @TokensSpec;

  $This->_GenerateLexer($Input, @TokensSpec);

  return $This;
}

# Generate a lexer using reference to an input iterator function, an open file
# handle or an input string passed as first parameter by the caller along
# with token specifications as second paramater...
#
sub _GenerateLexer {
  my($This, $Input, @TokensSpec) = @_;

  if ($This->{InputType} =~ /^InputIterator$/i) {
    $This->_GenerateInputIteratorLexer($Input, @TokensSpec);
  }
  elsif ($This->{InputType} =~ /^FileStream$/i) {
    $This->_GenerateInputFileStreamLexer($Input, @TokensSpec);
  }
  elsif ($This->{InputType} =~ /^String$/i) {
    $This->_GenerateInputStringLexer($Input, @TokensSpec);
  }
  else {
    croak "Error: ${ClassName}->new: Object can't be instantiated: Invalid input parameter type specified. Supported parameter types: a reference to input iterator function, a reference to an open file handle or a text string...";
  }

  return $This;
}

# Generate a lexer using specifed input iterator...
#
sub _GenerateInputIteratorLexer {
  my($This, $InputIteratorRef, @TokensSpec) = @_;

  $This->_GenerateChainedLexer($InputIteratorRef, @TokensSpec);

  return $This;
}

# Generate a lexer using specifed input file stream reference...
#
sub _GenerateInputFileStreamLexer {
  my($This, $FileHandleRef, @TokensSpec) = @_;

  # Iterator is a annoymous function reference and Perl keeps $FileHandleRef
  # in scope during its execution.

  $This->_GenerateChainedLexer( sub { return <$FileHandleRef>; }, @TokensSpec);

  return $This;
}

# Generate a lexer using specifed input string...
#
sub _GenerateInputStringLexer {
  my($This, $Text, @TokensSpec) = @_;
  my(@InputText) = ($Text);

  # Iterator is a annoymous function reference and Perl keeps @InputText
  # in scope during its execution.

  $This->_GenerateChainedLexer( sub { return shift @InputText; }, @TokensSpec);

  return $This;
}

# Get next available token label and value pair as an array reference or unrecognized
# text from input stream by either removing it from the input or simply peeking ahead...
#
# Supported mode values: Peek, Next. Default: Next
#
sub Lex {
  my($This, $Mode) = @_;

  return $This->{ChainedLexer}->($Mode)
}

# Get next available token label and value pair as an array reference or unrecognized
# text from input stream by either removing it from the input stream...
#
sub Next {
  my($This) = @_;

  return $This->Lex();
}

# Get next available token label and value pair as an array reference or unrecognized
# text from input stream by simply peeking ahead and without removing it from the input
# stream..
#
sub Peek {
  my($This) = @_;

  return $This->Lex('Peek')
}

# Get a reference to lexer method to be used by the caller...
#
sub GetLex {
  my($This) = @_;

  return sub { $This->Lex(); };
}

# The chained lexer generation is implemented based on examples in Higher-order Perl
# [ Ref 126 ] book.
#
# Generate a lexer using specified input iterator and chaining it with other lexers generated
# for all token specifications. The lexer generated for first token specification uses input
# iterator to retrieve any available input text; the subsequent chained lexeres for rest
# of the tokens use lexers generated for previous token specifications to get next input
# which might be unmatched input text or a reference to an array containing token and
# matched text pair.
#
sub _GenerateChainedLexer {
  my($This, $InputIteratorRef, @TokensSpec) = @_;
  my($TokenSpecRef, $ChainedLexer);

  $ChainedLexer = undef;
  for $TokenSpecRef (@TokensSpec) {
    $ChainedLexer = defined $ChainedLexer ? $This->_GenerateLexerForToken($ChainedLexer, @{$TokenSpecRef}) : $This->_GenerateLexerForToken($InputIteratorRef, @{$TokenSpecRef});
  }

  $This->{ChainedLexer} = $ChainedLexer;

  return $This;
}


# Generate a lexer using specifed token specification using specified input or
# input retrieved using another token lexer. The lexer retrieving input from the
# specified input stream is at the bottom of the chain.
#
sub _GenerateLexerForToken {
  my($This, $InputIteratorOrLexer, $TokenLabel, $RegexPattern, $TokenMatchActionRef) = @_;
  my($TokenMatchAndSplitRef, $InputBuffer, @ProcessedTokens);

  # Input buffer for a specific lexer in chained lexers containing unprocessed
  # text for token specifications retrieved from a downstrean lexer or intial
  # input...
  #
  $InputBuffer = "";

  # @ProcessedTokens contains either references to an array containing token label
  # and matched text or any unmatched input text string...
  #
  @ProcessedTokens = ();

  # Setup a default annoymous function reference to generate an array reference
  # containing $Token and text matched to $RegexPattern.
  #
  $TokenMatchActionRef = defined $TokenMatchActionRef ? $TokenMatchActionRef : sub { my($This, $Label, $MatchedText) = @_; return [$Label, $MatchedText]  };

  # Setup an annoymous function to match and split input text using $RegexPattern for
  # a specific token during its lexer invocation in chained lexers.
  #
  # The usage of parenthesis around $RegexPattern during split allows capturing of matched
  # text, which is subsequently processed to retrieve matched $Token values. The split function
  # inserts a "" separator in the returned array as first entry whenever $InputText starts with
  # $RegexPattern. $InputText is returned as the only element for no match.
  #
  $TokenMatchAndSplitRef = sub { my($InputText) = @_; return split /($RegexPattern)/, $InputText; };

  # Setup a lexer for $TokenLabel as an annoymous function and return its reference to caller
  # which in turns chains the lexers for all $Tokens before returning a reference to a lexer
  # at top of the lexer chain.
  #
  # Perl maintains scope of all variables defined with in the scope of the current function
  # during invocation of annoymous function even after the return call.
  #
  return sub {
    my($Mode) = @_;

    # Currenly supported value for mode: Peek, Next
    #
    $Mode = defined $Mode ? $Mode : 'Next';

    while (@ProcessedTokens == 0 && defined $InputBuffer ) {
      # Get any new input....
      my $NewInput = $InputIteratorOrLexer->();

      if (ref $NewInput) {
	# Input is an array reference containing matched token and text returned by
	# a chained lexer downstream lexer...
	#
	# Match $RegexPattern in available buffer text to retieve any matched text
	# for current $Token. $Separator might be "": $RegexPattern is at start of
	# of $InputBuffer
	#
	# Process input buffer containing text to be matched for the current lexer
	# which didn't get processed earlier during @NewTokens > 2  while loop:
	# no match for current lexer or more input available. It maintains order
	# of token matching in input stream.
	#
	my($Separator, $MatchedTokenRefOrText);

	($Separator, $MatchedTokenRefOrText) = $TokenMatchAndSplitRef->($InputBuffer);
	if (defined $MatchedTokenRefOrText) {
	  $MatchedTokenRefOrText = $TokenMatchActionRef->($This, $TokenLabel, $MatchedTokenRefOrText);
	}

	# Collect valid token references or text...
	push @ProcessedTokens, grep { defined $_ && $_ ne "" } ($Separator, $MatchedTokenRefOrText, $NewInput);

	# Empty put buffer...
	$InputBuffer = "";

	# Get out of the loop as processed token refererences and/or text  are available...
	last;
      }

      # Process input retrieved from downstream lexer or input iterator which hasn't
      # been processed into tokens..
      if (defined $NewInput) {
	$InputBuffer .= $NewInput;
      }

      # Retrieve any matched tokens from available input for the current lexer...
      #
      my(@NewTokens) = $TokenMatchAndSplitRef->($InputBuffer);

      while ( @NewTokens > 2 || @NewTokens && !defined $NewInput) {
	# Scenario 1: Complete match
	#   @NewTokens > 2 : Availability of separator, matched token text, separator.
	#   The separator might correspond to token for a token for upstream lexer followed
	#   by matched token from current lexer. It ends up getting passed to upsrteam
	#   lexer for processing.
	#
	# Scenario 2: No more input available from iterator or downstream lexer
	#   @NewTokens <= 2 and no more input implies any left over text in buffer. And
	#   it ends up getting passed to upsrteam for processing.
	#

	# Take off any unprocessed input text that doesn't match off the buffer: It'll be
	# passed to upstream chained lexer for processing...
	#
	push @ProcessedTokens, shift @NewTokens;

	if (@NewTokens) {
	  my $MatchedTokenText = shift @NewTokens;
	  push @ProcessedTokens, $TokenMatchActionRef->($This, $TokenLabel, $MatchedTokenText);
	}
      }

      # Retrieve any leftover text from NewTokens and put it back into InputBuffer for
      # processing by current lexer. All token references have been taken out....
      #
      $InputBuffer = "";
      if (@NewTokens) {
	$InputBuffer = join "", @NewTokens;
      }

      if (!defined $NewInput) {
	# No more input from the downstream lexer...
	$InputBuffer = undef;
      }

      # Clean up any empty strings from ProcessedTokens containing token
      # references or text...
      @ProcessedTokens = grep { $_ ne "" } @ProcessedTokens;

    }

    # Return reference to an array containing token and matched text or just unmatched input text...
    my $TokenRefOrText = undef;

    if (@ProcessedTokens) {
      # Get first available reference either by just peeking or removing it from the list
      # of available tokens...
      $TokenRefOrText = ($Mode =~ /^Peek$/i) ?  $ProcessedTokens[0] : shift @ProcessedTokens;
    }

    return $TokenRefOrText;
  };
}

# Is it a lexer object?
sub _IsLexer {
  my($Object) = @_;

  return (Scalar::Util::blessed($Object) && $Object->isa($ClassName)) ? 1 : 0;
}

# Return a string containing information about lexer...
sub StringifyLexer {
  my($This) = @_;
  my($LexerString);

  $LexerString = "Lexer: PackageName: $ClassName; " . $This->_GetLexerInfoString();

  return $LexerString;
}

# Return a string containing information about lexer...
sub _GetLexerInfoString {
  my($This) = @_;
  my($LexerInfoString, $TokensSpec, $TokenSpec, $TokenLabel, $TokenMatchRegex, $TokenMatchAction);

  $LexerInfoString = "InputType: $This->{InputType}";

  if ($This->{InputType} =~ /^String$/i) {
    $LexerInfoString .= "; InputString: $This->{Input}";
  }

  $TokensSpec = "TokensSpecifications: <None>";
  if (@{$This->{TokensSpec}}) {
    $TokensSpec = "TokensSpecifications: < [Label, MatchRegex, MatchAction]:";
    for $TokenSpec (@{$This->{TokensSpec}}) {
      ($TokenLabel, $TokenMatchRegex) = @{$TokenSpec};
      $TokenMatchAction = (@{$TokenSpec} == 3) ? "$TokenSpec->[2]" : "undefined";
      $TokensSpec .= " [$TokenLabel, $TokenMatchRegex, $TokenMatchAction]";
    }
    $TokensSpec .= " >";
  }

  $LexerInfoString .= "; $TokensSpec";

  return $LexerInfoString;
}

1;

__END__

=head1 NAME

Parsers::Lexer

=head1 SYNOPSIS

use Parsers::Lexer;

use Parsers::Lexer qw(:all);

=head1 DESCRIPTION

B<Lexer> class provides the following methods:

new, GetLex, Lex, Next, Peek, StringifyLexer

The object oriented chained B<Lexer> is implemented based on examples available in
Higher-order Perl [ Ref 126 ] book by Mark J. Dominus. It is designed to be used
both in standalone mode or as a base class for B<YYLexer>.

A chained lexer is created by generating a lexer for for the first specified token
specification using specified input and chaining it with other lexers generated for all
subsequent token specifications. The lexer generated for the first token specification
uses input iterator to retrieve any available input text; the subsequent chained lexeres
for rest of the token specifications use lexers generated for previous token
specifications to get next input, which might be unmatched input text or a reference
to an array containing token and  matched text pair.

=head2 METHODS

=over 4

=item B<new>

    $Lexer = new Parsers::Lexer($Input, @TokensSpec);

Using specified I<Input> and I<TokensSpec>, B<new> method generates a new lexer
and returns a reference to newly created B<Lexer> object.

Example:

    # Tokens specifications supplied by the caller. It's an array containing references
    # to arrays with each containing TokenLabel and TokenMatchRegex pair along with
    # an option reference to code to be executed after a matched.
    #
    @LexerTokensSpec = (
        [ 'LETTER', qr/[a-zA-Z]/ ],
        [ 'NUMBER', qr/\d+/ ],
        [ 'SPACE', qr/[ ]*/,
            sub { my($This, $TokenLabel, $MatchedText) = @_; return ''; }
        ],
        [ 'NEWLINE', qr/(?:\r\n|\r|\n)/,
            sub { my($This, $TokenLabel, $MatchedText) = @_;  return "\n"; }
        ],
        [ 'CHAR', qr/./ ]
    );

    # Input string...
    $InputText = 'y = 3 + 4';
    $Lexer = new Parsers::Lexer($InputText, @LexerTokensSpec);

    # Process input stream...
    while (defined($Token = $Lexer->Lex())) {
        print "Token: " . ((ref $Token) ? "@{$Token}" : "$Token") . "\n";
    }

    # Input file...
    $InputFile = "Input.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $Lexer = new Parsers::Lexer(\*INPUTFILE, @LexerTokensSpec);

    # Input file iterator...
    $InputFile = "TestSimpleCalcParser.txt";
    open INPUTFILE, "$InputFile" or die "Couldn't open $InputFile: $!\n";
    $InputIterator = sub { return <INPUTFILE>; };
    $Lexer = new Parsers::Lexer($InputIterator, @LexerTokensSpec);

    @LexerTokensSpec = (
        [ 'VAR', qr/[[:alpha:]]+/ ],
        [ 'NUM', qr/\d+/ ],
        [ 'OP', qr/[-+=\/]/,
            sub { my($This, $Label, $Value) = @_;
                $Value .= "; ord: " . ord $Value;
                return [$Label, $Value];
            }
        ],
        [ 'NEWLINE', qr/(?:\r\n|\r|\n)/, sub { return [$_[1], 'NewLine']; } ],
        [ 'SPACE', qr/\s*/, sub { return [$_[1], 'Space']; } ],
    );

    # Look ahead without removing...
    $Token = $Lexer->Lex('Peek');
    if (defined $Token && ref $Token) {
        print "PEEK: Token: @{$Token}\n\n";
    }

    # Process input stream...
    while (defined($Token = $Lexer->Lex())) {
        print "Token: " . ((ref $Token) ? "@{$Token}" : "$Token") . "\n";
    }

=item B<GetLex>

    $LexerRef = $Lexer->GetLex();

Returns a refernece to I<Lexer> method to the caller for use in a specific B<YYLexer>.

=item B<Lex>

    $TokenRefOrText = $Lexer->Lex($Mode);
    if (ref $TokenRefOrText) {
        ($TokenLabel, $TokenValue) = @{$TokenRefOrText};
    }
    else {
        $TokenText = $TokenRefOrText;
    }

Get next available token label and value pair as an array reference or unrecognized
text from input stream by either removing it from the input or simply peeking ahead
and without removing it from the input stream.

Possible I<Mode> values: I<Peek, Next>. Default: I<Next>.

=item B<Next>

    $TokenRefOrText = $Lexer->Next();

Get next available token label and value pair as an array reference or unrecognized
text from input stream by removing it from the input stream.

=item B<Peek>

    $TokenRefOrText = $Lexer->Peek();

Get next available token label and value pair as an array reference or unrecognized
text from input stream by by simply peeking ahead and without removing it from the
input stream.

=item B<StringifyLexer>

    $LexerString = $Lexer->StringifyLexer();

Returns a string containing information about I<Lexer> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

YYLexer.pm, SimpleCalcYYLexer.pm, SimpleCalcParser.yy

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
