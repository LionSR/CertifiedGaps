# Contributing

Anyone can add a result to the leaderboard of proven gaps. A result can be a proven gap for a new
model, a larger proven gap for a known point, or a Lean theorem for a bound that so far has a
computer proof. The [contribute page](https://sirui-lu.com/haldane-gap-diagrams/contribute.html)
of the website explains each case (source: [`home_page/contribute.md`](home_page/contribute.md)).
The [leaderboard](https://sirui-lu.com/haldane-gap-diagrams/leaderboard.html) lists the results
and the open challenges.

## In short

1. A computer proof goes in `models/<name>/` and follows [`AGENTS.md`](AGENTS.md). It needs an
   exact reflection-positivity identity, rigorous bounds at the starting scale, the gap bound
   from the repository's recurrence, and a checker program `check.jl`. Add a row to the summary
   table of the root [`README.md`](README.md).
2. A Lean theorem goes in the library `Leaderboard`, as described in
   [`leaderboard/README.md`](leaderboard/README.md). It may use only Lean's standard axioms.
3. Every result has a row in [`home_page/_data/records.yml`](home_page/_data/records.yml), whose
   header documents the fields. Check the file with

   ```sh
   julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
   julia --project=tools/leaderboard tools/leaderboard/check_records.jl
   ```

4. To claim an open challenge, use the "Take this on" button on the leaderboard. To propose a
   new model, use the "Propose a target model" issue form. To report an error, use the "Report a
   problem with a record" form.

## Checks

Automatic checks run on every pull request. They build the Lean library, list the axioms of every
theorem and check every proof again in the Lean kernel. They compare every leaderboard row with
its certificate and its Lean theorem in exact arithmetic, and they run the checker program of
each changed model. A reviewer then checks that each row describes the model its certificate and
its Lean theorem are about, and that its badge is right.
