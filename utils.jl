using WAV

const audacityExe = "/Applications/Audacity.app/Contents/MacOS/Audacity"

function audacity(x; fs=framerate(x), timeout=10.0)
  wavfname = tempname() * ".wav"
  wavwrite(samples(x), wavfname; Fs=fs)
  run(`"$audacityExe" "$wavfname"`; wait=false)
  @async begin
    sleep(timeout)
    rm(wavfname)
  end
  nothing
end

wavsize(filename) = wavread(filename; format="size")
