use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/lib", "$FindBin::Bin/..";

use Test::More;

use_ok('HGVSpOneLetter');

my $plugin = HGVSpOneLetter->new({});

is_deeply($plugin->feature_types, ['Transcript'], 'runs on transcripts only');
is_deeply([sort keys %{ $plugin->get_header_info }], ['HGVSp_OneLetter', 'HGVSp_Valid'], 'declares HGVSp_OneLetter and HGVSp_Valid');

# three-letter HGVSp as VEP 115 writes it => expected one-letter, expected HGVSp_Valid
# covers every HGVSp pattern seen in the simulated TP53/PTEN/GPX1 runs
my @cases = (
  # substitutions
  ['NP_000537.3:p.Arg175His',                   'NP_000537.3:p.R175H',              'true',  'missense'],
  ['NP_000305.3:p.Arg130Ter',                   'NP_000305.3:p.R130*',              'true',  'nonsense'],
  ['NP_000305.3:p.Asp173%3D',                   'NP_000305.3:p.D173%3D',            'true',  'synonymous, escaped'],
  ['NP_000537.3:p.Ter394%3D',                   'NP_000537.3:p.*394%3D',            'true',  'stop retained, escaped'],
  ['NP_000305.3:p.Met1?',                       'NP_000305.3:p.M1?',                'true',  'start lost'],
  # frameshifts
  ['NP_000305.3:p.Asp173GlyfsTer11',            'NP_000305.3:p.D173Gfs*11',         'true',  'frameshift'],
  ['NP_000572.2:p.Lys97GlnfsTer?',              'NP_000572.2:p.K97Qfs*?',           'true',  'frameshift, no new stop'],
  ['NP_000305.3:p.Lys6Ter',                     'NP_000305.3:p.K6*',                'true',  'frameshift with immediate stop (SingleLetterAA: *fs*1)'],
  ['NP_000537.3:p.Ter394CysfsTer77',            'NP_000537.3:p.*394Cfs*77',         'false', 'frameshift at the stop codon (HGVS: extension)'],
  ['NP_000572.2:p.Ter204TrpfsTer?',             'NP_000572.2:p.*204Wfs*?',          'false', 'frameshift at the stop codon, no new stop'],
  # deletions, duplications, insertions
  ['NP_000305.3:p.Met174del',                   'NP_000305.3:p.M174del',            'true',  'single deletion'],
  ['NP_000305.3:p.Met174_Thr175del',            'NP_000305.3:p.M174_T175del',       'true',  'multi deletion'],
  ['NP_000537.3:p.Asp393dup',                   'NP_000537.3:p.D393dup',            'true',  'single dup (SingleLetterAA: stays 3-letter)'],
  ['NP_000305.3:p.Pro172_Asp173dup',            'NP_000305.3:p.P172_D173dup',       'true',  'multi dup'],
  ['NP_000305.3:p.Pro172_Asp173insValPro',      'NP_000305.3:p.P172_D173insVP',     'true',  'insertion'],
  ['NP_000305.3:p.Met1_Thr2insTer',             'NP_000305.3:p.M1_T2ins*',          'true',  'insertion of stop'],
  ['NP_000537.3:p.Asp393_Ter394insCys',         'NP_000537.3:p.D393_*394insC',      'true',  'insertion before the stop'],
  ['NP_000305.3:p.Asp173delinsAlaHis',          'NP_000305.3:p.D173delinsAH',       'true',  'single delins'],
  ['NP_000305.3:p.Asp173_Met174delinsVal',      'NP_000305.3:p.D173_M174delinsV',   'true',  'multi delins'],
  ['NP_000305.3:p.Tyr27delinsTer',              'NP_000305.3:p.Y27delins*',         'true',  'delins to stop (SingleLetterAA: delins*D)'],
  # stop loss / extensions
  ['NP_000305.3:p.Ter404ArgextTer8',            'NP_000305.3:p.*404Rext*8',         'true',  'stop lost'],
  ['NP_000572.2:p.Ter204TyrextTer?',            'NP_000572.2:p.*204Yext*?',         'true',  'stop lost, no new stop'],
  ['NP_000305.3:p.Ter404delextTer7',            'NP_000305.3:p.*404delext*7',       'false', 'stop deleted (VEP form)'],
  ['NP_000305.3:p.Val403_Ter404delextTer6',     'NP_000305.3:p.V403_*404delext*6',  'false', 'stop deleted, multi (VEP form)'],
  ['NP_000305.3:p.Ter404delinsLeuArg',          'NP_000305.3:p.*404delinsLR',       'false', 'delins replacing stop, no ext (VEP form)'],
  ['NP_000305.3:p.Ter404delinsTer',             'NP_000305.3:p.*404delins*',        'false', 'stop replaced by stop (VEP form)'],
  # other VEP forms that are not valid HGVS: still converted code-for-code
  ['NP_000305.3:p.Met1_?2',                     'NP_000305.3:p.M1_?2',              'false', 'start-codon deletion (VEP form)'],
  ['NP_000305.3:p.MetThr1_?2',                  'NP_000305.3:p.MT1_?2',             'false', 'start-codon deletion, two residues (VEP form)'],
  ['NP_000305.3:p.LeuSer112%3D',                'NP_000305.3:p.LS112%3D',           'false', 'synonymous MNV across codons (VEP form)'],
  ['NP_000305.3:p.Val0_Val1insAsn',             'NP_000305.3:p.V0_V1insN',          'false', 'insertion before the start codon (VEP form)'],
  # rare residue codes
  ['NP_000572.2:p.Sec49Arg',                    'NP_000572.2:p.U49R',               'true',  'selenocysteine'],
  ['NP_000572.2:p.Sec49AlafsTer10',             'NP_000572.2:p.U49Afs*10',          'true',  'selenocysteine frameshift'],
  ['NP_000000.1:p.Pyl10Ala',                    'NP_000000.1:p.O10A',               'true',  'pyrrolysine'],
  ['NP_000000.1:p.Xaa10Ala',                    'NP_000000.1:p.X10A',               'true',  'unknown residue'],
  # flags that change VEP's HGVSp format
  ['NP_000537.3:p.(Arg175His)',                 'NP_000537.3:p.(R175H)',            'true',  '--hgvsp_use_prediction'],
  ['NP_000305.3:p.(Asp173=)',                   'NP_000305.3:p.(D173=)',            'true',  '--hgvsp_use_prediction --no_escape'],
  ['NP_000305.3:p.(Met1_?2)',                   'NP_000305.3:p.(M1_?2)',            'false', 'invalid form inside brackets'],
  ['NP_000537:p.Arg175His',                     'NP_000537:p.R175H',                'true',  '--remove_hgvsp_version'],
  ['ENSP00000269305.4:p.Arg175His',             'ENSP00000269305.4:p.R175H',        'true',  'Ensembl protein ID'],
);

for my $c (@cases) {
  my ($in, $want, $valid, $label) = @$c;
  is_deeply($plugin->run(undef, { HGVSp => $in }),
            { HGVSp_OneLetter => $want, HGVSp_Valid => $valid }, $label);
}

# every residue code converts and validates
my %aa = (
  Ala => 'A', Arg => 'R', Asn => 'N', Asp => 'D', Cys => 'C', Gln => 'Q', Glu => 'E',
  Gly => 'G', His => 'H', Ile => 'I', Leu => 'L', Lys => 'K', Met => 'M', Phe => 'F',
  Pro => 'P', Ser => 'S', Thr => 'T', Trp => 'W', Tyr => 'Y', Val => 'V',
  Sec => 'U', Pyl => 'O', Xaa => 'X', Ter => '*',
);
for my $three (sort keys %aa) {
  is_deeply($plugin->run(undef, { HGVSp => "NP_1.1:p.Ala5$three" }),
            { HGVSp_OneLetter => "NP_1.1:p.A5$aa{$three}", HGVSp_Valid => 'true' }, "$three -> $aa{$three}");
}

# syntax checks that VEP output has not (yet) exercised
my @invalid = (
  [':p.Arg175His',                  'missing protein ID'],
  ['NP_1.1:Arg175His',              'missing p.'],
  ['NP_1.1:p.',                     'empty description'],
  ['NP_1.1:p.Asx175His',            'unrecognised residue code (left three-letter)'],
  ['NP_1.1:p.Arg0His',              'position 0'],
  ['NP_1.1:p.Arg175',               'substitution with no new residue'],
  ['NP_1.1:p.Ter175His',            'stop substitution without ext'],
  ['NP_1.1:p.Arg175_Gly173del',     'range end before start'],
  ['NP_1.1:p.Arg175_Gly177insAla',  'insertion between non-adjacent residues'],
  ['NP_1.1:p.Arg175_Gly176ins',     'insertion with no sequence'],
  ['NP_1.1:p.Arg175delinsTerAla',   'residue after a stop in delins'],
  ['NP_1.1:p.Arg175_Gly176insTerAla','residue after a stop in insertion'],
  ['NP_1.1:p.Arg175TerfsTer3',      'frameshift to an immediate stop (should be substitution)'],
  ['NP_1.1:p.Arg175GlyfsTer0',      'frameshift with new stop at 0'],
  ['NP_1.1:p.Arg175GlyextTer5',     'extension not starting at a stop'],
  ['NP_1.1:p.(Arg175His',           'unbalanced bracket'],
  ['NP_1.1:p.Arg175His extra',      'trailing text'],
);
for my $c (@invalid) {
  my ($in, $label) = @$c;
  is($plugin->run(undef, { HGVSp => $in })->{HGVSp_Valid}, 'false', "invalid: $label");
}

my @valid = (
  ['NP_1.1:p.0',                    'no protein'],
  ['NP_1.1:p.?',                    'unknown effect'],
  ['NP_1.1:p.Arg175fs',             'short frameshift'],
  ['NP_1.1:p.Met1ext-5',            'N-terminal extension'],
  ['NP_1.1:p.(Arg175GlyfsTer?)',    'predicted frameshift, no new stop'],
);
for my $c (@valid) {
  my ($in, $label) = @$c;
  is($plugin->run(undef, { HGVSp => $in })->{HGVSp_Valid}, 'true', "valid: $label");
}

# nothing to convert
is_deeply($plugin->run(undef, {}), {}, 'no HGVSp (non-coding / no --hgvs)');
is_deeply($plugin->run(undef, { HGVSp => undef }), {}, 'undefined HGVSp');
is_deeply($plugin->run(undef, { HGVSp => 'p.Arg175His' }), {}, 'no protein ID separator');

# VEP's own HGVSp is not modified
my $line = { HGVSp => 'NP_000537.3:p.Arg175His' };
$plugin->run(undef, $line);
is($line->{HGVSp}, 'NP_000537.3:p.Arg175His', 'input HGVSp left unchanged');

done_testing();
