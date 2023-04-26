// prepare footer

// switch footer visibility
let footer = document.getElementsByTagName('footer')[0];
let footer_down = true;
let show_footer = document.getElementById("show_footer");
show_footer.addEventListener('click',(ev)=>{
  if (footer_down) {
    footer.style.transform = "translate(0)";
    show_footer.style.transform = "rotate(0.5turn)";
    footer_down = false;
  } else {
    footer.style.transform = "translate( 0, calc(100% - 100px) )";
    show_footer.style.transform = "rotate(0turn)";
    footer_down = true;
  }
});

// set double slider for order
var slider = document.getElementById('order_slider');
noUiSlider.create(slider, {
    start: [1, 6],
    step: 1,
    connect: true,
    tooltips: true,
    format: {
      to:   (value)=>{ return Math.ceil(Number(value)); },
      from: (value)=>{ return Math.ceil(Number(value)); }
    },
    range: {
        'min': 1,
        'max': 10
    }
});

// sorting and filter algorithms

// sort methods
const sort_methods = {
  'alpha': (elm) => { return elm.id; } ,
  'order': (elm) => { return Number(elm.getAttribute("data-order")); } ,
  'stages': (elm) => { return Number(elm.getAttribute("data-nstages")); } ,
  'stage_order': (elm) => { return Number(elm.getAttribute("data-stage_order")); } ,
  'xmax': (elm) => { return Number(elm.getAttribute("data-xmax")); } ,
  'ymax': (elm) => { return Number(elm.getAttribute("data-ymax")); } ,
  'random': (elm) => { return Math.random() - 0.5; }
}

let sort_btn = document.getElementById("sort_btn");
let sort_by_elm = document.getElementById("sort");

function sort_rklist( method ) {
  let frag = document.createDocumentFragment();
  let viewer = sort_methods[method];
  let sortedList = [...ul.children].sort((a,b) => {
    return  viewer(a) > viewer(b);
  });
  for ( const item of sortedList ) {
    frag.appendChild(item);
  }
  ul.appendChild(frag);
}

sort_btn.addEventListener('click',function (event) {
  sort_rklist(sort_by_elm.value);
});

// filter methods
const filter_methods = {
  'explicit': (elm) => { return elm.getAttribute("data-explicit") === "true";  } ,
  'implicit': (elm) => { return elm.getAttribute("data-explicit") === "false"; } ,
  'dirk':     (elm) => { return elm.getAttribute("data-dirk")     === "true";  } ,
  'embedded': (elm) => { return elm.getAttribute("data-embedded") === "true";  } ,
  'reset':    (elm) => { return true; } ,
  'order':    (elm) => { return elm.getAttribute("data-order"); }
};
let order_slider = document.getElementById("order_slider").noUiSlider;

let filter_btn = document.getElementById("filter_btn");
filter_btn.addEventListener('click',function (event) {
  let frag = document.createDocumentFragment();

  let filter_by_elm;
  let types = document.querySelectorAll("input[name='filter_type']");
  for ( const radio of types ) {
    if (radio.checked) {
      filter_by_elm = radio;
      break;
    }
  }

  const range_order = order_slider.get();

  let viewer = filter_methods[filter_by_elm.value];
  let filtredList = [...ul.children].map( (elm) => {
    elm.classList.remove('unselected');
    let o = filter_methods['order'](elm);
    if ( !( viewer(elm) && ( range_order[0] <= o && o <= range_order[1] ) ) ) {
      elm.classList.add('unselected');
    }
    return elm;
  });
  for ( const item of filtredList ) {
    frag.appendChild(item);
  }
  ul.appendChild(frag);
});
