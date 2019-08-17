function surveys = get_surveys()
% Which strata are part of which survey...

s = 1;
surveys(s).name = 'CCAMLR 2019';
surveys(s).strata = ["AP", "SS", "SSI", "SOI", "Sand", "ESS"]; s = s + 1;

surveys(s).name = 'AMLR 2019';
surveys(s).strata = ["Joinville", "Elephant", "West", "Bransfield"]; s = s + 1;

surveys(s).name = 'Norway1';
surveys(s).strata = "SOF"; s = s + 1;

surveys(s).name = 'Norway1';
surveys(s).strata = "SOC"; s = s + 1;
