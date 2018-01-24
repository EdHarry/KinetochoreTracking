function [mode, f, fjac, user] = objFcn_wrapper(mode, m, n, ldfj, needfi, x, fjac, nstate, user)

[mode, f, fjac, user] = objFcn_mex(int32(mode), int32(m), int32(n), int32(ldfj), int32(needfi), x, fjac, int32(nstate), user);

end

