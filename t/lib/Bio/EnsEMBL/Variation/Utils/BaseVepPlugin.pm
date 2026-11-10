package Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

# Minimal stand-in for the ensembl-variation base class, so the plugin can be
# unit tested without a VEP installation. Only what HGVSpOneLetter uses.

use strict;
use warnings;

sub new {
  my ($class, $config, @params) = @_;
  return bless { config => $config, params => \@params }, $class;
}

1;
