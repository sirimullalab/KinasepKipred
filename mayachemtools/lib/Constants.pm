package Constants;
#
# File: Constants.pm
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
use Exporter;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);

# Groups of constants...
my(@MathConstants) = qw(Pi TwoPi PiByTwo PiByFour InversePi OneByPi SqrtTwoPi OneBySqrtTwoPi InverseSqrtTwoPi Phi E BigNumber);

# Export all constants...
@EXPORT = (@MathConstants);
@EXPORT_OK = qw();

# Tags to export individual groups of constants...
%EXPORT_TAGS = (
		math => [@MathConstants],
		all  => [@EXPORT, @EXPORT_OK]
	       );

# Math constants..
use constant {
  Pi => CORE::atan2(0, -1),
  TwoPi => 2 * CORE::atan2(0, -1),
  PiByTwo => CORE::atan2(1, 0),
  PiByFour => CORE::atan2(1, 1),
  OneByPi => 1.0 / CORE::atan2(0, -1),
  InversePi => 1.0 / CORE::atan2(0, -1),
  SqrtTwoPi => CORE::sqrt(2 * CORE::atan2(0, -1)),
  OneBySqrtTwoPi => 1 / CORE::sqrt(2 * CORE::atan2(0, -1)),
  InverseSqrtTwoPi => 1 / CORE::sqrt(2 * CORE::atan2(0, -1)),

  Phi => (1 + CORE::sqrt(5))/2, # Golden ratio...

  E => CORE::exp(1),

  BigNumber => 2**31 - 1, # Unsigned 31 bit number

};

1;

__END__

=head1 NAME

Constants

=head1 SYNOPSIS

use Constants;

use Constants qw(:all);

use Constants qw(:math);

=head1 DESCRIPTION

B<Constants> module provided the following constants:

Pi, TwoPi, PiByTwo, PiByFour, InversePi, OneByPi, SqrtTwoPi, OneBySqrtTwoPi, InverseSqrtTwoPi,
Phi, E

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

MathUtil.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
