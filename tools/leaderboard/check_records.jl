# Validate the leaderboard data home_page/_data/records.yml (schema in its header comment):
#
#     julia --project=tools/leaderboard tools/leaderboard/check_records.jl
#     julia --project=tools/leaderboard tools/leaderboard/check_records.jl \
#         --lean-records home_page/_data/lean_records.json      # also compare the Lean records
#
# Exits nonzero on any error. See `LeaderboardCheck.main` for the options.

using LeaderboardCheck

exit(LeaderboardCheck.main(ARGS))
