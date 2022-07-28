function add_btn(list_id) {
  let btn_add = document.createElement("button");
  btn_add.id = "add";
  btn_add.innerHTML = "<span class='add'>+</span>";
  btn_add.addEventListener('click',()=>{

    let meth = document.createElement("article");
    meth.classList.add('method');
    
    let h2 = label("rk_"+document.getElementsByClassName('method').length,"New method");
    h2.setAttribute('contenteditable',true);
    meth.appendChild(h2);

    // first choose number of stages
    let input = document.createElement("input");
    input.type = 'number';
    input.placeholder = "Number of stages";
    meth.appendChild(input);

    let btn = document.createElement("button");
    btn.appendChild(document.createTextNode("Make tableau!"));
    btn.classList.add("maketable-btn");

    // when click "make tableau!" change element
    btn.addEventListener('click',(ev)=>{
      let this_btn = ev.target;
      let this_parent = this_btn.parentNode;
      let nstages_input = this_parent.children[1];

      const Nstages = Number( nstages_input.value );
      if ( Nstages > 0 ) {
        // create editable empty Butcher tableau
        let tab = document.createElement("table");
        tab.classList.add('butcher');
        for(i=0;i<=Nstages;++i) {
          let tr = document.createElement("tr");
          for(j=0;j<=Nstages;++j) {
            let td = document.createElement("td");
            td.setAttribute('contenteditable',true);
            if (i!=Nstages) {
              if ( j == 0 ) {
                td.innerHTML = "c_" + (i+1);
              } else {
                td.innerHTML = "a_{" + (i+1) + "," + j + "}";
              }
            } else {
              if ( j == 0 ) {
                td.innerHTML = "";
                td.setAttribute('contenteditable',false);
              } else {
                td.innerHTML = "b_" + j;
              }
            }
            tr.appendChild(td);
          }
          tab.appendChild(tr);
        }
        this_parent.appendChild(tab);

        let action_btn = document.createElement("button");
        action_btn.appendChild(document.createTextNode("Serialize..."));
        action_btn.classList.add("serialize-btn");
        action_btn.addEventListener("click",(ev)=>{
          let this_btn = ev.target;
          let this_parent = this_btn.parentNode;

          let label = this_parent.children[0].innerText.trim();
          let A = Array(Nstages);
          let b = Array(Nstages);
          let c = Array(Nstages);
          
          // read Butcher tableau (and parse expressions)
          let tab = this_parent.children[1];
          for (i=0;i<tab.childElementCount-1;++i) {
            c[i] = tab.children[i].children[0].innerText.trim().replaceAll("^","**");
            A[i] = Array(Nstages);
            for (j=1;j<tab.childElementCount;++j) {
              A[i][j-1] = tab.children[i].children[j].innerText.trim().replaceAll("^","**");
            }
          }
          for (i=1;i<tab.childElementCount;++i) {
            b[i-1] = tab.children[tab.childElementCount-1].children[i].innerText.trim().replaceAll("^","**");
          }

          let data = {'label':label,'A':A,'b':b,'c':c};

          let isopen = false;
          // if already a details section, remove it ! (then rebuild it)
          if ( this_parent.lastChild.localName == "details" ) {
            isopen = this_parent.lastChild.open;
            this_parent.removeChild(this_parent.lastChild);
          }
          let visu = document.createElement("details");
          visu.classList.add("visu");
          visu.open = isopen;
          let save_zone = document.createElement("textarea");
          save_zone.value = JSON.stringify(data,null,2);
          visu.appendChild(save_zone);
          
          let preview = document.createElement("div");
          preview.classList.add("preview");
          visu.appendChild(rk_to_elm(data,preview));

          this_parent.appendChild(visu);
        });
        this_parent.appendChild(action_btn);

        // if N > 0 remove tableau generator
        this_parent.removeChild(nstages_input);
        this_parent.removeChild(this_btn);
      }
    });
    meth.appendChild(btn);
    
    document.getElementById(list_id).appendChild(meth);
    meth.scrollIntoView();
  });

  return btn_add;
}

document.getElementById("btn_add").appendChild(add_btn("list"));
