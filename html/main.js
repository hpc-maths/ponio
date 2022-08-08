const parser = new AsciiMathParser(); // to convert python expression (asciimath ich) into LaTeX expression
function asciimath(expr) {
  return expr.replaceAll("**","^").replaceAll(/\d\*\d/ig,"x");
}
function render(expr,elm) {
  katex.render(parser.parse(asciimath(expr)),elm,{displayMode:true});
  return elm;
}

function label(id,lab) {
  let header = document.createElement("header");
  //let h2 = document.createElement("h2");
  //h2.appendChild(document.createTextNode(lab));
  //header.appendChild(h2);
  let a = document.createElement("a");
  a.href = "#" + id;
  a.appendChild(document.createTextNode(lab));
  header.appendChild(a);
  
  return header;
}

function butcher_tableau(A,b,c,b2=undefined) {
  let is_embedded = ( b2 !== undefined );

  let container = document.createElement("div");
  container.classList.add("butcher_container");
  let tab = document.createElement("table");
  tab.classList.add("butcher");
  if (is_embedded) { tab.classList.add("embedded"); }
  for (i=0;i<A.length;++i) {
    let tr = document.createElement("tr");
    let ci = document.createElement("td");
    render(c[i],ci);
    tr.appendChild(ci);
    for (j=0;j<A[i].length;++j) {
      let aij = document.createElement("td");
      render(A[i][j],aij);
      tr.appendChild(aij);
    }
    tab.appendChild(tr);
  }

  let tr = document.createElement("tr");
  tr.appendChild(document.createElement("td"));
  for (i=0;i<b.length;++i) {
    let bi = document.createElement("td");
    render(b[i],bi);
    tr.appendChild(bi);
  }
  tab.appendChild(tr);

  if ( is_embedded ) {
    let tr = document.createElement("tr");
    tr.appendChild(document.createElement("td"));
    for (i=0;i<b2.length;++i) {
      let b2i = document.createElement("td");
      render(b2[i],b2i);
      tr.appendChild(b2i);
    }
    tab.appendChild(tr);
  }

  container.appendChild(tab);

  return container;
}

function stability_function(R_expr) {
  let div = document.createElement("div");
  div.classList.add("stability_function");
  render(R_expr,div);
  
  return div;
}
function resume_tableau(nstages,order,stage_order) {
  let tab = document.createElement("table");
  tab.classList.add("resume_tableau");
  
  let tr1 = document.createElement("tr");
  let th1 = document.createElement("th"); th1.setAttribute('scope',"row"); th1.appendChild(document.createTextNode("#stages"));
  let td1 = document.createElement("td"); td1.appendChild(document.createTextNode(nstages));
  tr1.appendChild(th1); tr1.appendChild(td1);
  tab.appendChild(tr1);

  let tr2 = document.createElement("tr");
  let th2 = document.createElement("th"); th2.setAttribute('scope',"row"); th2.appendChild(document.createTextNode("order"));
  let td2 = document.createElement("td"); td2.appendChild(document.createTextNode(order));
  tr2.appendChild(th2); tr2.appendChild(td2);
  tab.appendChild(tr2);

  let tr3 = document.createElement("tr");
  let th3 = document.createElement("th"); th3.setAttribute('scope',"row"); th3.appendChild(document.createTextNode("stage order"));
  let td3 = document.createElement("td"); td3.appendChild(document.createTextNode(stage_order));
  tr3.appendChild(th3); tr3.appendChild(td3);
  tab.appendChild(tr3);

  return tab;
}

function expr_to_canvas(f,[xmin,xmax],[ymin,ymax],thresholds) {
  let canvas = document.createElement("canvas");
  canvas.height = 300;
  canvas.width  = 300;
  const n = 512, m = 512;

  context = canvas.getContext("2d");
  projection = d3.geoIdentity().scale(canvas.width / n);
  path = d3.geoPath(projection,context);

  let values = new Array(n * m);
  for (let j = 0, k = 0; j < m; ++j) {
    for (let i = 0; i < n; ++i, ++k) {
      values[k] = f( i * (xmax-xmin)/ n + xmin, j * (ymax-ymin)/ m + ymin);
    }
  }

  const contours = d3.contours()
    .size([n,m])
    .thresholds(thresholds)
    (values);
  contours.forEach( c => {
    if (c.coordinates.length == 0) return;
    context.beginPath();
    path(c);
    context.fillStyle = "#fefefe";
    context.fill();
    //context.stroke();
  });

  return canvas;
}
function data_to_canvas(data,[xmin,xmax],[ymin,ymax],thresholds,options={},outside=true,data2=undefined) {
  let canvas = document.createElement("canvas");
  canvas.height = 512;
  canvas.width  = 512;
  const n = data[0].length , m = data.length;

  context = canvas.getContext("2d");
  projection = d3.geoIdentity().scale(canvas.width / n);
  path = d3.geoPath(projection,context);

  const contours = d3.contours()
    .size([n,m])
    .thresholds(thresholds)
    ( [].concat(...data) );
  contours.forEach( c => {
    if (c.coordinates.length == 0) return;
    context.beginPath();
    path(c);
    if ( outside ) {
      context.rect(0,0,canvas.width,canvas.height);
    }
    context.fillStyle = options.hasOwnProperty('fill_color') ? options.fill_color : "#ff0000";
    context.fill();
  });
  if ( data2 !== undefined ) {
    const contours2 = d3.contours()
      .size([n,m])
      .thresholds(thresholds)
      ( [].concat(...data2) );
    contours2.forEach( c => {
      if (c.coordinates.length == 0) return;
      context.beginPath();
      path(c);
      context.lineWidth = 3.0;
      context.strokeStyle = options.hasOwnProperty('stroke_color') ? options.stroke_color : "#000000";
      context.stroke();
    });
  }

  // draw axis
  const dx = (xmax-xmin)/canvas.width;
  const dy = (ymax-ymin)/canvas.height;
  const i_0 = parseInt((-xmin)/dx), j_0 = parseInt(-ymin/dy);

  context.strokeStyle = options.hasOwnProperty('axis_color') ? options.axis_color : "#000000";
  context.lineWidth = 1.0;
  // real axis
  context.beginPath();
  context.moveTo(0,j_0);
  context.lineTo(n,j_0);
  context.closePath();
  context.stroke();
  // imaginary axis
  context.beginPath();
  context.moveTo(i_0,0);
  context.lineTo(i_0,n);
  context.closePath();
  context.stroke();

  // xticks
  for ( x=parseInt(xmin)-1 ; x<xmax ; ++x ) {
    context.beginPath();
    context.moveTo(parseInt((x-xmin)/dx),j_0-5);
    context.lineTo(parseInt((x-xmin)/dx),j_0+5);
    context.closePath();
    context.stroke();
  }
  // yticks
  for ( y=parseInt(ymin)-1 ; y<ymax ; ++y ) {
    context.beginPath();
    context.moveTo(i_0-5,parseInt((y-ymin)/dy));
    context.lineTo(i_0+5,parseInt((y-ymin)/dy));
    context.closePath();
    context.stroke();
  }

  
  return canvas;
}

function stability_domain( data , options ) {
  let fig = document.createElement("figure");
  fig.appendChild( data_to_canvas(data.data,[data.xmin,data.xmax],[data.ymin,data.ymax],[1.0],options=options,outside=true,data2=data.data_embedded) );
  let caption = document.createElement("figcaption");
  let def = document.createElement("span"); katex.render(String.raw`\mathcal{D} = \{ z\in\mathbb{C}\,/\,|R(z)| \leq 1 \}`,def,{displayMode:true});
  caption.innerHTML = "Stability domain: " + def.outerHTML;
  fig.appendChild(caption);
  return fig;
}
function order_star( data , options ) {
  let fig = document.createElement("figure");
  fig.appendChild( data_to_canvas(data.data,[data.xmin,data.xmax],[data.ymin,data.ymax],[1.0],options=options) );
  let caption = document.createElement("figcaption");
  let def = document.createElement("span"); katex.render(String.raw`\mathcal{A}_- = \{ z\in\mathbb{C}\,/\,|e^{-z}R(z)| < 1 \}`,def,{displayMode:true});
  caption.innerHTML = "Order star: " + def.outerHTML;
  fig.appendChild(caption);
  return fig;
}

function prepare_canvas(id) {
  let div = document.createElement("div");
  div.classList.add("drawing");
  div.setAttribute("data-id",id);

  return div;
}

function scheme(stages) {
  let ul = document.createElement("ul");
  ul.classList.add("scheme");
  stages.map(eq => {
    let li = document.createElement("li");
    katex.render(eq,li,{displayMode:true});
    ul.appendChild(li);
  });
  return ul;
}
function code(c,id) {
  let pre = document.createElement("pre");
  let code = document.createElement("code");
  code.classList.add("language-py");
  code.innerHTML = hljs.highlight(c,{language: 'python'}).value;
  code.id = "code-"+id;
  
  let btn = document.createElement('button');
  btn.classList.add("clipboard");
  btn.setAttribute("data-clipboard-target","#"+code.id);
  btn.innerHTML = "<svg xmlns='http://www.w3.org/2000/svg' width='20' height='20' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='2' stroke-linecap='round' stroke-linejoin='round' class='feather feather-paperclip'><path d='M21.44 11.05l-9.19 9.19a6 6 0 0 1-8.49-8.49l9.19-9.19a4 4 0 0 1 5.66 5.66l-9.2 9.19a2 2 0 0 1-2.83-2.83l8.49-8.48'></path></svg><span>Copy code</span>";

  pre.appendChild(btn);
  pre.appendChild(code);

  return pre;
}
function lscheme(stages) {
  let ul = document.createElement("ul");
  ul.classList.add("lawson-scheme");
  stages.map(eq => {
    let li = document.createElement("li");
    katex.render(eq,li,{displayMode:true});
    ul.appendChild(li);
  });
  return ul;
}

function rk_to_elm(rk,elm,options) {
  elm.id = rk.id;
  elm.classList.add('method');
  elm.setAttribute("data-order",rk.order);
  elm.setAttribute("data-nstages",rk.nstages);
  elm.setAttribute("data-stage_order",rk.stage_order);
  elm.setAttribute("data-xmax",0.);
  elm.setAttribute("data-ymax",0.);
  elm.setAttribute("data-explicit",rk.is_explicit);
  elm.setAttribute("data-dirk",rk.is_dirk);
  elm.setAttribute("data-embedded",rk.is_embedded);
  
  elm.appendChild(label(rk.id,rk.label));
  elm.appendChild(butcher_tableau(rk.A,rk.b,rk.c,rk.b2));
  
  if (rk.hasOwnProperty('stability_function')) {
    let details = document.createElement("details");
    
    let summary = document.createElement("summary");
    summary.appendChild(document.createTextNode("Details"))
    details.appendChild(summary);
    
    let p = document.createElement("p");
    let title_R = document.createElement("h6");
    let R = document.createElement("span"); katex.render("R(z)",R);
    title_R.innerHTML = "Stability function "+R.outerHTML+":";
    p.appendChild(title_R);
    p.appendChild(stability_function(rk.stability_function));
    p.appendChild(resume_tableau(rk.nstages,rk.order,rk.stage_order));
    p.appendChild(prepare_canvas(rk.id));
    
    summary.addEventListener('click',function (event){
        let div = event.target.parentElement.children[1].getElementsByClassName("drawing")[0];
        let id = div.getAttribute("data-id");
        const rk_url = new URL(window.location.href.split("/").slice(0,-1).join("/") + "/api/" + rk.id + ".json");
        fetch(rk_url)
          .then( res => res.json() )
          .then( out => {
            div.appendChild(stability_domain(out.stability_domain, options));
            div.appendChild(order_star(out.order_star, options));
          });
      },{'once':true});

    

    let title_edo = document.createElement("h6"); title_edo.appendChild(document.createTextNode("Runge-Kutta scheme:")); p.appendChild(title_edo);
    let edo = document.createElement("p"); katex.render("\\dot{u}=f(t,u)",edo,{displayMode:true});
    p.appendChild(edo);
    p.appendChild(scheme(rk.scheme));
  
    if (rk.hasOwnProperty('code')) {
      p.appendChild(code(rk.code,rk.id));
    }
    
    let title_ledo = document.createElement("h6"); title_ledo.appendChild(document.createTextNode("Lawson scheme:")); p.appendChild(title_ledo);
    let ledo = document.createElement("p"); katex.render("\\dot{u}=Lu + N(t,u)",ledo,{displayMode:true});
    p.appendChild(ledo);
    p.appendChild(scheme(rk.lawson_scheme));
  
    if (rk.hasOwnProperty('lawson_code')) {
      p.appendChild(code(rk.lawson_code,"l_"+rk.id));
      let warning = document.createElement("blockquote");
      warning.classList.add("warning");
      warning.innerHTML = "default NumPy function <a href='https://numpy.org/doc/stable/reference/generated/numpy.exp.html'><samp>np.exp</samp></a> for exponantial take the exponential of each terms of <samp>ndarray</samp>, which is not the exponential of a matrix needed by Lawson method. If this is for matrix computation, prefer <a href='https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html'><samp>scipy.linalg.expm</samp></a> function.";
      p.appendChild(warning);
    }

    if (rk.hasOwnProperty('doi')) {
      let doi = document.createElement("blockquote");
      doi.classList.add("doi");

      
      let url_crossref = "https://api.crossref.org/works/" + rk.doi;
      fetch(url_crossref)
      .then( res => res.json() )
      .then( res => res.message )
      .then( out => {
          let link = document.createElement("a");
          link.href = out.URL;
          link.innerHTML = "<span>Link to reference</span>";
          doi.appendChild(link);

          let authors = document.createElement("ul");
          authors.classList.add("authors");
          out.author.forEach( (author) => {
            let auth = document.createElement("li");
            auth.appendChild(document.createTextNode(author.family + ", " + author.given));
            authors.appendChild(auth);
          });
          doi.appendChild(authors);
          
          let title = document.createElement("span");
          title.classList.add("title");
          title.appendChild(document.createTextNode(out.title[0]));
          doi.appendChild(title);

          let pubdate = document.createElement("span");
          pubdate.classList.add("date");
          pubdate.appendChild(document.createTextNode(out.published['date-parts'][0][0]));
          doi.appendChild(pubdate);

          let publisher = document.createElement("span");
          publisher.classList.add("publisher");
          publisher.appendChild(document.createTextNode(out['short-container-title'][0]));
          doi.appendChild(publisher);
        });

        p.appendChild(doi);
    }

    details.appendChild(p);
    elm.appendChild(details);
  }

  return elm;
}