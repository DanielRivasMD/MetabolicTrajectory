module FFTStage

export run_fft

using Logging

function run_fft(params)
  @info "[fft] computing FFT"
  sleep(0.5)
  @info "[fft] done"
end

end
