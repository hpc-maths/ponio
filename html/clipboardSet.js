// set clipboard
clipboard = new ClipboardJS('.clipboard');
clipboard.on('success', function(e) {
    const reset_text = e.trigger.innerText;
    e.trigger.children[1].innerText = "Copied!";
    setTimeout( ()=>{ e.trigger.children[1].innerText = reset_text; } , 1000 );

    e.clearSelection();
});
