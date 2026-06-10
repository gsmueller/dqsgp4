% report_check(name, tf) — print PASS/FAIL for a symbolic check.
function report_check(name, tf)
  if tf
    fprintf('  [PASS]  %s\n', name);
  else
    fprintf('  [FAIL]  %s\n', name);
  end
end
