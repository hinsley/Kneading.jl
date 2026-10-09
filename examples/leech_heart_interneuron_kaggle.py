import argparse
import base64
import io
import json
import os
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BUNDLE = ["Project.toml", "src", "ext", "examples/leech_heart_interneuron_kneading.jl"]

RUNNER = r'''
import base64, io, os, subprocess, tarfile
PAYLOAD = "@@PAYLOAD@@"
ENV = @@ENV@@
WORK = "/kaggle/working"
ROOT = "/kaggle/tmp/Kneading.jl"
PROJECT = "/kaggle/tmp/leech"
os.makedirs(ROOT, exist_ok=True)
os.makedirs(PROJECT, exist_ok=True)
tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD)), mode="r:gz").extractall(ROOT, filter="data")

def sh(cmd, **kw):
    print("+", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True, **kw)

if not os.path.exists("/kaggle/tmp/julia/bin/julia"):
    sh("cd /kaggle/tmp && curl -fsSL -o julia.tgz https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.9-linux-x86_64.tar.gz"
       " && tar xzf julia.tgz && mv julia-1.11.9 julia")
JL = "/kaggle/tmp/julia/bin/julia"
sh("nvidia-smi")
if os.path.exists(ROOT + "/environment/Manifest.toml"):
    for name in ("Project.toml", "Manifest.toml"):
        text = open(f"{ROOT}/environment/{name}").read().replace("@@KNEADING@@", ROOT)
        open(f"{PROJECT}/{name}", "w").write(text)
    sh(f"{JL} --project={PROJECT} -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using CUDA; CUDA.versioninfo()'")
else:
    sh(f"{JL} --project={PROJECT} -e 'using Pkg; Pkg.develop(path=\"{ROOT}\"); "
       "Pkg.add([\"CUDA\", \"KernelAbstractions\", \"DynamicalSystemsBase\", \"StaticArrays\", \"ForwardDiff\"]); "
       "Pkg.precompile(); using CUDA; CUDA.versioninfo()'")
env = dict(os.environ)
env.update(ENV)
env["LEECH_OUTPUT"] = WORK + "/leech-heart-interneuron"
script = ROOT + "/examples/leech_heart_interneuron_kneading.jl"
call = "main()" if env.get("LEECH_DEVICE") == "cpu" else (
    "focus, = main(backend = CUDABackend()); "
    "problems = [FlowKneadingProblem(leech(r.metadata.parameters...); initializer = r.initialization, settings...) "
    "for r in focus.results if !isnothing(r)]; "
    "println(\"Device words alone: \", @elapsed(flow_kneading(problems; backend = CUDABackend())), \" s\")")
with open(WORK + "/scan.log", "w") as log:
    code = subprocess.run([JL, "-t", "auto", f"--project={PROJECT}", "-e",
        f"using CUDA; include(\"{script}\"); {call}"],
        env=env, stdout=log, stderr=subprocess.STDOUT).returncode
print(open(WORK + "/scan.log").read()[-6000:], flush=True)
if code:
    raise SystemExit(f"scan failed with exit code {code}")
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("username")
    parser.add_argument("slug", nargs="?", default="leech-heart-interneuron-kneading")
    parser.add_argument("--env", nargs="*", default=[])
    parser.add_argument("--out", default=os.path.join(REPO, "output", "kaggle"))
    parser.add_argument("--environment")
    args = parser.parse_args()
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for name in BUNDLE:
            archive.add(os.path.join(REPO, name), arcname=name)
        if args.environment:
            for name in ("Project.toml", "Manifest.toml"):
                text = open(os.path.join(args.environment, name)).read().replace(REPO, "@@KNEADING@@").encode()
                info = tarfile.TarInfo(f"environment/{name}")
                info.size = len(text)
                archive.addfile(info, io.BytesIO(text))
    environment = dict(item.split("=", 1) for item in args.env)
    script = RUNNER.replace("@@PAYLOAD@@", base64.b64encode(buffer.getvalue()).decode())
    script = script.replace("@@ENV@@", repr(environment))
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "run_scan.py"), "w") as handle:
        handle.write(script)
    metadata = {
        "id": f"{args.username}/{args.slug}",
        "title": args.slug,
        "code_file": "run_scan.py",
        "language": "python",
        "kernel_type": "script",
        "is_private": True,
        "enable_gpu": True,
        "machine_shape": "NvidiaTeslaT4",
        "enable_internet": True,
        "dataset_sources": [],
        "competition_sources": [],
        "kernel_sources": [],
    }
    with open(os.path.join(args.out, "kernel-metadata.json"), "w") as handle:
        json.dump(metadata, handle, indent=2)
    print(f"Wrote {args.out}/run_scan.py and kernel-metadata.json for {metadata['id']}")


if __name__ == "__main__":
    main()
