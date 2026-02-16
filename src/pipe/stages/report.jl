module ReportStage

export run_report

using Logging

function run_report(params)
  @info "[report] generating report"
  sleep(0.5)
  @info "[report] done"
end

end
