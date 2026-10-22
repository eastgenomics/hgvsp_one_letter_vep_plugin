=head1 LICENSE

 Licensed under the Apache License, Version 2.0. Structure based on the
 Ensembl VEP plugin SingleLetterAA (Copyright EMBL-European Bioinformatics
 Institute, Apache License 2.0).

=head1 NAME

 HGVSpOneLetter

=head1 VERSION

 1.0.0

=head1 SYNOPSIS

 cp HGVSpOneLetter.pm ~/.vep/Plugins/
 ./vep -i variations.vcf --hgvs --plugin HGVSpOneLetter

=head1 DESCRIPTION

 VEP plugin that adds an HGVSp_OneLetter field: the HGVSp with single-letter
 amino acid codes (Ter becomes *). VEP's own HGVSp field is left unchanged.

 It converts the HGVSp string VEP has already built for the output line
 (three-letter, via the default path of hgvs_protein()) code by code, rather
 than calling hgvs_protein() in single-letter mode, which mis-formats some
 stop-gain, insertion and duplication cases in Ensembl release 115.

 Because it starts from VEP's own HGVSp, --hgvsp_use_prediction,
 --remove_hgvsp_version and --no_escape are honoured. Requires --hgvs.

 It also adds HGVSp_Valid ("true"/"false"): whether VEP's three-letter HGVSp
 is a well-formed HGVS protein description. This checks syntax only (protein
 ID, residue codes, positions, and the forms of substitution, synonymous,
 deletion, duplication, insertion, delins, frameshift, extension and start
 loss); it does not check the change against the protein sequence. Strings
 that are not valid HGVS, such as VEP's p.Met1_?2, p.Ter404delextTer7 or
 p.LeuSer112=, and any unrecognised residue code, are flagged "false".
 HGVSp_OneLetter is still written for them.

=cut

package HGVSpOneLetter;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

our $VERSION = '1.0.0';

my %AA = (
  Ala => 'A', Arg => 'R', Asn => 'N', Asp => 'D', Cys => 'C', Gln => 'Q', Glu => 'E',
  Gly => 'G', His => 'H', Ile => 'I', Leu => 'L', Lys => 'K', Met => 'M', Phe => 'F',
  Pro => 'P', Ser => 'S', Thr => 'T', Trp => 'W', Tyr => 'Y', Val => 'V',
  Sec => 'U', Pyl => 'O', Xaa => 'X', Ter => '*',
);
my $AA_RE = join '|', keys %AA;

# residues other than stop, and residues including stop
my $RES  = qr/(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Sec|Pyl|Xaa)/;
my $RESX = qr/(?:$RES|Ter)/;
my $POS  = qr/([1-9][0-9]*)/;
# new stop position for fs / ext: number or unknown
my $NEWSTOP = qr/Ter(?:[1-9][0-9]*|\?)/;

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    HGVSp_OneLetter => 'HGVSp with single-letter amino acid codes (converted from three-letter HGVSp)',
    HGVSp_Valid     => 'Whether HGVSp is a well-formed HGVS protein description (true/false; syntax only)',
  };
}

sub run {
  my ($self, $tva, $line) = @_;

  # VEP has already built HGVSp (three-letter) in the line hash
  my $hgvsp = $line->{HGVSp};
  return {} unless defined $hgvsp;

  # convert only the description after the protein ID
  my ($id, $desc) = split /:/, $hgvsp, 2;
  return {} unless defined $desc;
  my $valid = valid_hgvsp($id, $desc) ? 'true' : 'false';
  $desc =~ s/($AA_RE)/$AA{$1}/g;

  return { HGVSp_OneLetter => "$id:$desc", HGVSp_Valid => $valid };
}

# Syntax check of a three-letter HGVS protein description (after the ':')
sub valid_hgvsp {
  my ($id, $desc) = @_;

  return 0 unless length $id && $id =~ /^[A-Za-z0-9_.]+\z/;
  return 0 unless $desc =~ s/^p\.//;

  # predicted format p.(...) and VEP's URI escaping of '='
  $desc = $1 if $desc =~ /^\((.+)\)\z/;
  $desc =~ s/%3D/=/g;

  # no protein / unknown effect / start lost
  return 1 if $desc eq '0' || $desc eq '?' || $desc eq 'Met1?';

  # synonymous: one residue, including a retained stop
  return 1 if $desc =~ /^$RESX$POS=\z/;

  # substitution (missense or nonsense); a substituted stop must be an extension
  return 1 if $desc =~ /^$RES$POS$RESX\z/;

  # deletion / duplication of one residue or a range (stop deletion must be an extension)
  return 1 if $desc =~ /^$RES$POS(?:del|dup)\z/;
  return $2 > $1 if $desc =~ /^$RES${POS}_$RES$POS(?:del|dup)\z/;

  # insertion between adjacent residues; a stop can only end the inserted sequence
  return $2 == $1 + 1 if $desc =~ /^$RES${POS}_$RESX${POS}ins(?:$RES)*$RESX\z/;

  # deletion-insertion; a stop can only end the inserted sequence
  return 1 if $desc =~ /^$RES${POS}delins(?:$RES)*$RESX\z/;
  return $2 > $1 if $desc =~ /^$RES${POS}_$RES${POS}delins(?:$RES)*$RESX\z/;

  # frameshift: long form with new stop, or short form (a frameshift at the stop codon is an extension)
  return 1 if $desc =~ /^$RES$POS(?:${RES}fs$NEWSTOP|fs)\z/;

  # extension: C-terminal (stop lost) or N-terminal (new upstream start)
  return 1 if $desc =~ /^Ter$POS${RES}ext$NEWSTOP\z/;
  return 1 if $desc =~ /^Met1ext-[1-9][0-9]*\z/;

  return 0;
}

1;
