using CertifiedGaps.Results, JSON3
root=@__DIR__
r=JSON3.read(read(joinpath(root,"results","summary_D4_1.json"),String))
x=only(x for x in r.rows if x.n0==6 && x.beta==3)
data=Dict("D"=>"4//1","beta"=>"3//1","n0"=>6,"bound"=>x.bound,
    "quoted_survey_bound"=>0.543842,"difference"=>x.bound-0.543842,
    "certified"=>false,
    "interpretation"=>"The numerical rerun does not reproduce the quoted survey value to six decimals. The quoted value is not used as a certificate input; complete dense enclosures are computed independently.",
    "source_hashes"=>source_hashes(root;dirs=(),files=("model.jl","scan.jl","survey_check.jl")))
write_result(joinpath(root,"results","survey-crosscheck.json"),data)
println(data)
