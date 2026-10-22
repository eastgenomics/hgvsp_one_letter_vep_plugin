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

 It also adds HGVSp_Valid ("true"/"false"): whether the biocommons hgvs
 package (1.5.7) can parse VEP's three-letter HGVSp, using a port of its
 protein grammar. This is a syntax check only; it does not check positions,
 ranges or the change against the protein sequence. VEP forms that hgvs
 rejects, such as p.Met1_?2, p.Ter404delextTer7 and p.LeuSer112=, are flagged
 "false", as are residue codes hgvs does not support (Pyl, Asx, Glx, and Xaa
 at a position). HGVSp_OneLetter is still written for them.

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

# Protein rules of the biocommons hgvs 1.5.7 grammar (hgvs/_data/hgvs.pymeta),
# used for HGVSp_Valid. The grammar is a PEG: a choice takes the first
# alternative that matches and never goes back to try the others, and repeats
# are greedy. Choices are therefore atomic groups (?>...) and repeats are
# possessive, so this accepts exactly what hgvs parses. Asx and Glx are in
# the grammar but hgvs then fails converting them (KeyError), so they are left
# out here; Pyl is not in the grammar.
my $AA3     = qr/(?:Ala|Cys|Asp|Glu|Phe|Gly|His|Ile|Lys|Leu|Met|Asn|Pro|Gln|Arg|Ser|Thr|Val|Trp|Tyr|Xaa|Sec)/;
my $AA1     = qr/[ACDEFGHIKLMNPQRSTVWYBZXU]/;
my $AA13    = qr/(?>$AA3|$AA1)/;
my $TERM13  = qr/(?>Ter|[X*])/;
my $AAT13   = qr/(?>Ter|$AA3|[X*]|$AA1)/;                                  # aat3 | aat1
my $AAT13_SEQ = qr/(?>Ter|(?:$AA3)++(?:Ter)?+|[X*]|(?:$AA1)++[X*]?+)/;      # aat3_seq | aat1_seq
my $NUM     = qr/[0-9]++/;
my $OFFSET  = qr/(?>[0-9]++|\?)?+/;                                        # fsext_offset
my $FS      = qr/fs(?:$TERM13$OFFSET)?+/;
my $EXT     = qr/ext(?:(?>$TERM13$OFFSET|(?:$AA13)?+-$NUM))?+/;
# pro_edit: fs | ext | subst | delins | ins | del | dup | ident
my $EDIT    = qr/(?>(?:$AAT13)?+$FS|(?:$AAT13)?+$EXT|$AAT13|\?|delins$AAT13_SEQ|ins$AAT13_SEQ|del|dup|=)/;
my $POS     = qr/(?>$TERM13|$AA13)$NUM/;                                   # def_p_pos
my $DEF_IVL = qr/(?>${POS}_$POS|$POS)/;                                    # def_p_interval
my $IVL     = qr/(?>$DEF_IVL|\($DEF_IVL\))/;                               # p_interval
# p_posedit, including p_posedit_special
my $POSEDIT = qr/(?>$IVL$EDIT|\($IVL$EDIT\)|=|\(=\)|0\?|0|\?)/;
my $NAME    = qr/(?>[A-Za-z0-9]|[-_](?=[A-Za-z0-9]))/;
my $ACCN    = qr/[A-Za-z](?:$NAME++)?+(?:\.[0-9]++)?+/;
my $GENE    = qr/(?:\([A-Za-z]$NAME++\))?+/;                               # opt_gene_expr

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
    HGVSp_Valid     => 'Whether HGVSp parses as HGVS with biocommons hgvs 1.5.7 (true/false; syntax only)',
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
  my $valid = valid_hgvsp($hgvsp) ? 'true' : 'false';
  $desc =~ s/($AA_RE)/$AA{$1}/g;

  return { HGVSp_OneLetter => "$id:$desc", HGVSp_Valid => $valid };
}

# Whether biocommons hgvs 1.5.7 parses the HGVSp (syntax only)
sub valid_hgvsp {
  my ($hgvsp) = @_;

  # VEP's URI escaping of '='
  (my $s = $hgvsp) =~ s/%3D/=/g;

  return $s =~ /^$ACCN$GENE:p\.$POSEDIT\z/ ? 1 : 0;
}

1;
