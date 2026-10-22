use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/lib", "$FindBin::Bin/..";

use Test::More;

use_ok('HGVSpOneLetter');

my $plugin = HGVSpOneLetter->new({});

is_deeply($plugin->feature_types, ['Transcript'], 'runs on transcripts only');
is_deeply([sort keys %{ $plugin->get_header_info }], ['HGVSp_OneLetter', 'HGVSp_SyntaxFlag'], 'declares HGVSp_OneLetter and HGVSp_SyntaxFlag');

# three-letter HGVSp as VEP 115 writes it => expected one-letter, expected HGVSp_SyntaxFlag ('' = not set)
# covers every HGVSp pattern seen in the simulated TP53/PTEN/GPX1 runs
my @cases = (
  # substitutions
  ['NP_000537.3:p.Arg175His',                    'NP_000537.3:p.R175H',              '',                 'missense'],
  ['NP_000305.3:p.Arg130Ter',                    'NP_000305.3:p.R130*',              '',                 'nonsense'],
  ['NP_000305.3:p.Asp173%3D',                    'NP_000305.3:p.D173%3D',            '',                 'synonymous, escaped'],
  ['NP_000537.3:p.Ter394%3D',                    'NP_000537.3:p.*394%3D',            '',                 'stop retained, escaped'],
  ['NP_000305.3:p.Met1?',                        'NP_000305.3:p.M1?',                '',                 'start lost'],
  # frameshifts
  ['NP_000305.3:p.Asp173GlyfsTer11',             'NP_000305.3:p.D173Gfs*11',         '',                 'frameshift'],
  ['NP_000572.2:p.Lys97GlnfsTer?',               'NP_000572.2:p.K97Qfs*?',           '',                 'frameshift, no new stop'],
  ['NP_000305.3:p.Lys6Ter',                      'NP_000305.3:p.K6*',                '',                 'frameshift with immediate stop (SingleLetterAA: *fs*1)'],
  ['NP_000537.3:p.Ter394CysfsTer77',             'NP_000537.3:p.*394Cfs*77',         '',                 'frameshift at the stop codon (HGVS prefers ext; hgvs parses it)'],
  ['NP_000572.2:p.Ter204TrpfsTer?',              'NP_000572.2:p.*204Wfs*?',          '',                 'frameshift at the stop codon, no new stop'],
  # deletions, duplications, insertions
  ['NP_000305.3:p.Met174del',                    'NP_000305.3:p.M174del',            '',                 'single deletion'],
  ['NP_000305.3:p.Met174_Thr175del',             'NP_000305.3:p.M174_T175del',       '',                 'multi deletion'],
  ['NP_000537.3:p.Asp393dup',                    'NP_000537.3:p.D393dup',            '',                 'single dup (SingleLetterAA: stays 3-letter)'],
  ['NP_000305.3:p.Pro172_Asp173dup',             'NP_000305.3:p.P172_D173dup',       '',                 'multi dup'],
  ['NP_000305.3:p.Pro172_Asp173insValPro',       'NP_000305.3:p.P172_D173insVP',     '',                 'insertion'],
  ['NP_000305.3:p.Met1_Thr2insTer',              'NP_000305.3:p.M1_T2ins*',          '',                 'insertion of stop'],
  ['NP_000537.3:p.Asp393_Ter394insCys',          'NP_000537.3:p.D393_*394insC',      '',                 'insertion before the stop'],
  ['NP_000305.3:p.Asp173delinsAlaHis',           'NP_000305.3:p.D173delinsAH',       '',                 'single delins'],
  ['NP_000305.3:p.Asp173_Met174delinsVal',       'NP_000305.3:p.D173_M174delinsV',   '',                 'multi delins'],
  ['NP_000305.3:p.Tyr27delinsTer',               'NP_000305.3:p.Y27delins*',         '',                 'delins to stop (SingleLetterAA: delins*D)'],
  # stop loss / extensions
  ['NP_000305.3:p.Ter404ArgextTer8',             'NP_000305.3:p.*404Rext*8',         '',                 'stop lost'],
  ['NP_000572.2:p.Ter204TyrextTer?',             'NP_000572.2:p.*204Yext*?',         '',                 'stop lost, no new stop'],
  ['NP_000305.3:p.Ter404delextTer7',             'NP_000305.3:p.*404delext*7',       'hgvs_unparseable', 'stop deleted (VEP form)'],
  ['NP_000305.3:p.Val403_Ter404delextTer6',      'NP_000305.3:p.V403_*404delext*6',  'hgvs_unparseable', 'stop deleted, multi (VEP form)'],
  ['NP_000305.3:p.Ter404delinsLeuArg',           'NP_000305.3:p.*404delinsLR',       '',                 'delins replacing stop, no ext (VEP form; hgvs parses it)'],
  ['NP_000305.3:p.Ter404delinsTer',              'NP_000305.3:p.*404delins*',        '',                 'stop replaced by stop (VEP form; hgvs parses it)'],
  # VEP forms hgvs cannot parse: still converted code-for-code
  ['NP_000305.3:p.Met1_?2',                      'NP_000305.3:p.M1_?2',              'hgvs_unparseable', 'start-codon deletion (VEP form)'],
  ['NP_000305.3:p.MetThr1_?2',                   'NP_000305.3:p.MT1_?2',             'hgvs_unparseable', 'start-codon deletion, two residues (VEP form)'],
  ['NP_000305.3:p.LeuSer112%3D',                 'NP_000305.3:p.LS112%3D',           'hgvs_unparseable', 'synonymous MNV across codons (VEP form)'],
  ['NP_000305.3:p.Val0_Val1insAsn',              'NP_000305.3:p.V0_V1insN',          '',                 'insertion before the start codon (VEP form; hgvs parses it)'],
  # rare residue codes
  ['NP_000572.2:p.Sec49Arg',                     'NP_000572.2:p.U49R',               '',                 'selenocysteine'],
  ['NP_000572.2:p.Sec49AlafsTer10',              'NP_000572.2:p.U49Afs*10',          '',                 'selenocysteine frameshift'],
  ['NP_000000.1:p.Pyl10Ala',                     'NP_000000.1:p.O10A',               'hgvs_unparseable', 'pyrrolysine (not in the hgvs grammar)'],
  ['NP_000000.1:p.Xaa10Ala',                     'NP_000000.1:p.X10A',               'hgvs_unparseable', 'unknown residue at a position (hgvs reads X as a stop)'],
  # flags that change VEP's HGVSp format
  ['NP_000537.3:p.(Arg175His)',                  'NP_000537.3:p.(R175H)',            '',                 '--hgvsp_use_prediction'],
  ['NP_000305.3:p.(Asp173=)',                    'NP_000305.3:p.(D173=)',            '',                 '--hgvsp_use_prediction --no_escape'],
  ['NP_000305.3:p.(Met1_?2)',                    'NP_000305.3:p.(M1_?2)',            'hgvs_unparseable', 'invalid form inside brackets'],
  ['NP_000537:p.Arg175His',                      'NP_000537:p.R175H',                '',                 '--remove_hgvsp_version'],
  ['ENSP00000269305.4:p.Arg175His',              'ENSP00000269305.4:p.R175H',        '',                 'Ensembl protein ID'],
);

for my $c (@cases) {
  my ($in, $want, $flag, $label) = @$c;
  is_deeply($plugin->run(undef, { HGVSp => $in }),
            { HGVSp_OneLetter => $want, $flag ? (HGVSp_SyntaxFlag => $flag) : () }, $label);
}

# every residue code converts; all but Pyl are in the hgvs grammar
my %aa = (
  Ala => 'A', Arg => 'R', Asn => 'N', Asp => 'D', Cys => 'C', Gln => 'Q', Glu => 'E',
  Gly => 'G', His => 'H', Ile => 'I', Leu => 'L', Lys => 'K', Met => 'M', Phe => 'F',
  Pro => 'P', Ser => 'S', Thr => 'T', Trp => 'W', Tyr => 'Y', Val => 'V',
  Sec => 'U', Pyl => 'O', Xaa => 'X', Ter => '*',
);
for my $three (sort keys %aa) {
  is_deeply($plugin->run(undef, { HGVSp => "NP_1.1:p.Ala5$three" }),
            { HGVSp_OneLetter => "NP_1.1:p.A5$aa{$three}", $three eq 'Pyl' ? (HGVSp_SyntaxFlag => 'hgvs_unparseable') : () },
            "$three -> $aa{$three}");
}

# HGVSp_SyntaxFlag follows the biocommons hgvs 1.5.7 parser: these are strings it rejects
my @flagged = (
  [':p.Arg175His',                      'missing protein ID'],
  ['NP_1.1:Arg175His',                  'missing p.'],
  ['NP_1.1:c.175del',                   'not a protein description (hgvs parses c.; HGVSp must be p.)'],
  ['NP_1.1:p.',                         'empty description'],
  ['NP_1.1:p.Asx175His',                'ambiguous residue (hgvs KeyError)'],
  ['NP_1.1:p.Ala175Glx',                'ambiguous residue as the new residue'],
  ['NP_1.1:p.Arg175',                   'substitution with no new residue'],
  ['NP_1.1:p.Arg175_Gly176ins',         'insertion with no sequence'],
  ['NP_1.1:p.Arg175delinsTerAla',       'residue after a stop in delins'],
  ['NP_1.1:p.Arg175_Gly176insTerAla',   'residue after a stop in insertion'],
  ['NP_1.1:p.(Arg175His',               'unbalanced bracket'],
  ['NP_1.1:p.Arg175His extra',          'trailing text'],
  ["NP_1.1:p.Arg175His\n",              'trailing newline'],
  ['NP_1.:p.Arg175His',                 'accession version with no number'],
  ['1NP:p.Arg175His',                   'accession not starting with a letter'],
  ['NP_1.1(A):p.Arg175His',             'one-character gene symbol'],
);
for my $c (@flagged) {
  my ($in, $label) = @$c;
  is($plugin->run(undef, { HGVSp => $in })->{HGVSp_SyntaxFlag}, 'hgvs_unparseable', "flagged: $label");
}

# not checked by hgvs: positions, ordering and stop placement are syntax-only
my @not_flagged = (
  ['NP_1.1:p.0',                        'no protein'],
  ['NP_1.1:p.?',                        'unknown effect'],
  ['NP_1.1:p.Arg175fs',                 'short frameshift'],
  ['NP_1.1:p.Met1ext-5',                'N-terminal extension'],
  ['NP_1.1:p.(Arg175GlyfsTer?)',        'predicted frameshift, no new stop'],
  ['NP_1.1:p.Cys244?',                  'unknown effect mid-protein'],
  ['NP_1.1:p.Arg0His',                  'position 0'],
  ['NP_1.1:p.Ter175His',                'stop substitution without ext'],
  ['NP_1.1:p.Arg175_Gly173del',         'range end before start'],
  ['NP_1.1:p.Arg175_Gly177insAla',      'insertion between non-adjacent residues'],
  ['NP_1.1:p.Arg175TerfsTer3',          'frameshift to an immediate stop'],
  ['NP_1.1:p.Arg175GlyfsTer0',          'frameshift with new stop at 0'],
  ['NP_1.1:p.Arg175GlyextTer5',         'extension not starting at a stop'],
  ['NP_1.1(TP53):p.Arg175His',          'gene symbol after the accession'],
  ['NP-1:p.Arg175His',                  'hyphen in accession'],
  ['NP_1.1:p.=',                        'no change'],
  ['NP_1.1:p.(=)',                      'predicted no change'],
  ['NP_1.1:p.0?',                       'predicted no protein'],
  ['NP_1.1:p.R175*',                    'one-letter codes, * stop'],
  ['NP_1.1:p.R175X',                    'one-letter codes, X stop'],
  ['NP_1.1:p.A5insAV*',                 'one-letter insertion ending in a stop'],
  ['NP_1.1:p.A5delinsVX',               'one-letter delins ending in X'],
  ['NP_1.1:p.Arg175Glyfs*?',            'frameshift with one-letter stop'],
);
for my $c (@not_flagged) {
  my ($in, $label) = @$c;
  is($plugin->run(undef, { HGVSp => $in })->{HGVSp_SyntaxFlag}, undef, "not flagged: $label");
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
