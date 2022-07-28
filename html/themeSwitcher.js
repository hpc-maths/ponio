// change theme
let switcher = document.getElementById("switcher");
let light = true;
switcher.addEventListener('click',(ev)=>{
  html = document.getElementsByTagName("html")[0];
  if (light) {
    html.setAttribute("data-theme","dark")
    light = false;
    switcher.children[0].innerHTML = "🐈";
  } else {
    html.setAttribute("data-theme","light")
    light = true;
    switcher.children[0].innerHTML = "🐈‍⬛";
  }
});
