

function C = mmult(A,B)


ArBr = mmx_simple(real(A),real(B));
AiBi = mmx_simple(imag(A),imag(B));
ArpAiBrpBi = mmx_simple(real(A)+imag(A),real(B)+imag(B));

C = (ArBr - AiBi) + 1i*(ArpAiBrpBi - ArBr - AiBi);

end