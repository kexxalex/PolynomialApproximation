import subprocess

DIR = '/'.join(__file__.replace('\\', '/').split('/')[:-1]) # integration for jupyter notebooks

def approx(coeff: iter, data: list, digits=2, engine: str = "Fast.exe"):
    """
    coeff: iterable of coefficients to fit:
        (0, 2, 3) would fit the polynomial a3 * x^3 + a2 * x^2 + a_0
        (1,) would fit the polynomial a1 * x
    data: measured points of the form [(x,y), (x,y), (x,y), ...]
    digits: rounds the fit coefficients
    engine:
        "CLANG.exe" uses the clang++ compilation
        ".exe" the VisualStudio
    """
    args = [f"{DIR}/PolyApprox{engine}", str(len(coeff)), str(len(data)), *map(str, coeff)]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, creationflags=0x08000000)

    d = '\n'.join(map(lambda p: f"{float(p[0])}\n{float(p[1])}", data)).encode()
    out = proc.communicate(input=d)[0]

    if proc.returncode == 0:
        variance, *coefficients = out.decode("utf-8").split(';')[::-1]
        return {"variance": round(float(variance), digits*2), "coeff": tuple(map(lambda x: round(float(x), digits), coefficients))}

    elif proc.returncode == 1:
        print("PolyFitter needs more arguments")
    elif proc.returncode == 2:
        print("Coefficient count is 0")
    elif proc.returncode == 3:
        print("PolyFitter needs more fitting data")
    return {"variance": 0.0, "coeff": (0,)}
