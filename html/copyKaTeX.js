// simple way to copy equation in LaTeX mode
document.addEventListener('copy', function (event) {
  const selection = window.getSelection();
  if (selection.isCollapsed) {
    return;  // default action OK
  }
  const fragment = selection.getRangeAt(0).cloneContents();
  const katexs = fragment.querySelectorAll('.katex');
  if (katexs.length === 0) {
    return;  // default action OK
  }
  katexs.forEach (function (element) {
    const texSource = element.querySelector('annotation');
    if (texSource) {
      element.replaceWith(texSource);
      texSource.innerHTML = '$' + texSource.innerHTML + '$';
    }
  });
  fragment.querySelectorAll('.katex-display annotation').forEach (function (element) {
    element.innerHTML = '$' + element.innerHTML + '$';
  })
  event.clipboardData.setData('text/plain', fragment.textContent);
  event.clipboardData.setData('text/html', selection.getRangeAt(0).cloneContents());
  event.preventDefault();
})
